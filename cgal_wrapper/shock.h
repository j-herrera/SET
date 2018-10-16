// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#include <algorithm>
#include <math.h>
#include <random>
#include "types.h"
#include "./PolyhedralSurf.h"
#include "econe.h"
#include "path.h"
#include "interp.h"
#include "curvature.h"
#include "detachment.h"
#include "dbscan.h"

#include <CGAL/intersections_d.h>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include <CGAL/bounding_box.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>

#include <CGAL/convex_hull_3.h>

#include <CGAL/Side_of_triangle_mesh.h>

typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef CGAL::cpp11::array<std::size_t,3> Face;

template <class HDS>
class BuildPoly : public CGAL::Modifier_base<HDS> {

const std::vector<Point> *pts;
const std::vector<Face> *fcs;

public:
BuildPoly(const std::vector<Point> *_pts, const std::vector<Face> *_fcs) : pts(_pts), fcs(_fcs) {
}

void operator() ( HDS& hds) {
  CGAL::Polyhedron_incremental_builder_3<HDS> B(hds);
  B.begin_surface( (*pts).size(), (*fcs).size());

  for(std::vector<Point>::const_iterator it = (*pts).begin(); it != (*pts).end(); ++it ) {
    B.add_vertex(*it);
  }

  for( std::vector<Face>::const_iterator it = (*fcs).begin(); it != (*fcs).end(); ++it ) {
    if(B.test_facet((*it).begin(), (*it).end())) {
      B.add_facet((*it).begin(), (*it).end());
    } else {
      std::array< std::size_t, 3 > nf { {(*it)[0], (*it)[2], (*it)[1]} };
      if(B.test_facet(nf.begin(), nf.end())) {
        B.add_facet(nf.begin(), nf.end());
      } else {
        // std::cerr << "Error during facet insertion" << std::endl;
      }
    }
  }

  if(B.error())
  {
    std::cerr << "An error occured while creating a Polyhedron" << std::endl;
    B.rollback();
  }

  B.end_surface();
}
};

bool pcomp (Point i, Point j) {
  return ((i-j).squared_length() < 1.e-4);
}

std::vector<Polyhedron> compute_shock(Polyhedron & G, double gamma, double Mach, Vector & u, int model, double cut_angle, bool multi_shock, double eps, int minpts) {
  // Noise generator
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-1.e-12,1.e-12);

  // Compute streamlines and stagnation points
  Vector mu = -u;
  output op = path(G, mu, cut_angle);
  std::vector<std::vector<Point> > pstreams = op.streams;

  std::vector<std::vector<Point> > streams;
  std::map<Point,std::pair<Point, Point> > start_map;

  // Streamline breaking on compression
  for (unsigned int k = 0; k < pstreams.size(); k++) {
    std::vector<Point> stline = pstreams[k];
    std::vector<std::vector<Point> > substlines;

    Vector pface(stline.rbegin()[1]-stline.rbegin()[0]);
    pface = pface / sqrt(pface.squared_length());

    double prev_angle = acos(mu * pface);
    long unsigned int pc = stline.size();
    if(stline.size() > 3) {
      for (unsigned int i = stline.size()-2; i > 0; i--) {
        Vector nface(stline[i]-stline[i+1]);
        nface = nface / sqrt(nface.squared_length());
        double cangle = acos(mu * nface);
        if ((cangle * 0.7 > prev_angle) && cangle > 1.e-3) {
          prev_angle = cangle;
          std::vector<Point> tmpline;
          tmpline.insert (tmpline.begin(), &stline[i], &stline[std::min(pc+1,stline.size())]);
          substlines.push_back(tmpline);
          unsigned int lasti = i;

          for (unsigned int j = i-1; j > 0; j--) {
            nface = (stline[j]-stline[j+1]);
            nface = nface / sqrt(nface.squared_length());
            cangle = acos(mu * nface);
            i = j;

            if (cangle <= prev_angle) {
              start_map.emplace(stline[j], std::make_pair(stline[lasti], stline[lasti+1]));
              prev_angle = cangle;
              pc = j;
              break;
            } else {
              prev_angle = cangle;
            }
          }
        } else {
          prev_angle = cangle;
        }
      }
    }

    std::vector<Point> tmpline;
    tmpline.insert (tmpline.begin(), &stline[0], &stline[std::min(pc+1,stline.size())]);
    substlines.push_back(tmpline);

    for (unsigned int l = 0; l < substlines.size(); l++) {
      std::vector<Point> stlinet = substlines[l];
      if (stlinet.size()<3) {
        continue;
      }
      streams.push_back(stlinet);
    }
  }

  std::list<Point> st_points;
  for(unsigned int i = 0; i < streams.size(); i++) {
    st_points.push_back(streams[i].back());
    std::reverse(streams[i].begin(),streams[i].end());
  }
  st_points.unique();

  // Compute geometry limits
  Kernel::Iso_cuboid_3 c3 = CGAL::bounding_box(G.points_begin(), G.points_end());
  double d3 = sqrt((c3.max() - c3.min()).squared_length());
  std::pair<double, double> maxAng = maxAngle(model, Mach, gamma);
  double ma = maxAng.first;

  int nedg = 0;
  double tlen = 0.0;
  for ( Facet_iterator i = G.facets_begin(); i != G.facets_end(); ++i ) {
    Halfedge_handle h = i->halfedge();
    std::vector<double> edgs;
    edgs.push_back((h->vertex()->point() - h->opposite()->vertex()->point()).squared_length());
    edgs.push_back((h->next()->vertex()->point() - h->next()->opposite()->vertex()->point()).squared_length());
    edgs.push_back((h->prev()->vertex()->point() - h->prev()->opposite()->vertex()->point()).squared_length());

    tlen += sqrt(*std::min_element(edgs.begin(), edgs.end()));
    nedg++;
  }

  tlen /= nedg;

  // Compute stagnation clusters
  int clusters = 1;
  std::vector<std::pair<Point, int> > D;
  if (multi_shock) {
    std::vector<Point> st_pt;
    for(unsigned int i = 0; i < streams.size(); i++) {
      st_pt.push_back(streams[i][0]);
    }

    Kernel::Iso_cuboid_3 c_st_pt = CGAL::bounding_box(st_pt.begin(), st_pt.end());

    for(unsigned int i = 0; i < streams.size(); i++) {
      double ptx = (streams[i][0].x() - (c_st_pt.xmax()+c_st_pt.xmin())/2.0) / (c_st_pt.xmax() - c_st_pt.xmin());
      double pty = (streams[i][0].y() - (c_st_pt.ymax()+c_st_pt.ymin())/2.0) / (c_st_pt.ymax() - c_st_pt.ymin());
      double ptz = (streams[i][0].z() - (c_st_pt.zmax()+c_st_pt.zmin())/2.0) / (c_st_pt.zmax() - c_st_pt.zmin());
      std::pair<Point, int> ptm(Point(ptx, pty, ptz), -2);
      D.push_back(ptm);
    }
    clusters = DBSCAN(D, eps, minpts);
  }

  // Order clusters
  std::vector<int> cluster_order;
  if (multi_shock) {
    Vector iv = u;
    Vector jv = CGAL::cross_product( iv, Vector(std::rand(), std::rand(), std::rand()));
    jv = jv / sqrt(jv.squared_length());
    Vector kv = CGAL::cross_product( iv, jv);

    auto comp = [&](const std::pair<Point, int> &a, std::pair<Point, int> &b) -> bool {
                  if ( (a.first-Point(0.0, 0.0, 0.0)) * iv != (b.first-Point(0.0, 0.0, 0.0)) * iv) {
                    return (a.first-Point(0.0, 0.0, 0.0)) * iv > (b.first-Point(0.0, 0.0, 0.0)) * iv;
                  }

                  if ( (a.first-Point(0.0, 0.0, 0.0)) * jv != (b.first-Point(0.0, 0.0, 0.0)) * jv) {
                    return (a.first-Point(0.0, 0.0, 0.0)) * jv > (b.first-Point(0.0, 0.0, 0.0)) * jv;
                  }

                  return (a.first-Point(0.0, 0.0, 0.0)) * kv > (b.first-Point(0.0, 0.0, 0.0)) * kv;
                };

    std::vector<std::pair<Point, int> > DD = D;
    std::sort (DD.begin(), DD.end(), comp);

    for(unsigned int i = 0; i < DD.size(); i++) {
      if (DD[i].second == -1) {
        continue;
      }
      bool flag = true;
      for (unsigned int j=0; j< cluster_order.size(); j++) {
        if(DD[i].second == cluster_order[j]) {
          flag = false;
          break;
        }
      }
      if (flag) {
        cluster_order.push_back(DD[i].second);
      }
    }
  } else{
    cluster_order.push_back(1);
  }

  // Compute curvatures at stagnation points
  std::map< Point, std::pair<std::array<double, 2>, std::vector<Vector> > > ks;
  ks = curvatures(G, st_points);

  // Compute curvature map
  std::list<Point> vtx;
  for ( Point_iterator v = G.points_begin(); v != G.points_end(); ++v ) {
    vtx.push_back(*v);
  }

  std::map< Point, std::pair<std::array<double, 2>, std::vector<Vector> > > kmap;
  kmap = curvatures(G, vtx);

  // Compute AABB tree
  Tree tree(faces(G).first, faces(G).second, G);
  tree.accelerate_distance_queries();

  // Compute point-normal pairs over streamlines
  std::vector<Polyhedron> S;
  std::vector<Polyhedron> Sch;
  std::vector<double> after_mach;
  for (std::vector<int>::iterator k = cluster_order.begin(); k < cluster_order.end(); k++) {
    std::vector<Point> points;
    double local_mach = Mach;
    for (unsigned int i = 0; i < streams.size(); i++) {
      if (multi_shock && *k != D[i].second) {
        continue;
      }

      Point start = streams[i][0];

      if (multi_shock) {
        local_mach = Mach;
        std::vector<Polyhedron>::reverse_iterator rit = S.rbegin();
        std::vector<Polyhedron>::reverse_iterator rkt = Sch.rbegin();
        std::vector<double>::reverse_iterator rjt = after_mach.rbegin();
        for (; rit!= S.rend(); ++rit, ++rjt, ++rkt) {
          CGAL::Side_of_triangle_mesh<Polyhedron, Kernel> inside(*rkt);

          CGAL::Bounded_side res = inside(streams[i][0]);

          if (res == CGAL::ON_BOUNDED_SIDE || res == CGAL::ON_BOUNDARY) {
            local_mach = *rjt;
            break;
          }
        }
        //recompute ma

        try {
          std::pair<Point, Point> ppair = start_map.at(streams[i][0]);
          Vector prface = ppair.first-ppair.second;
          prface = prface / sqrt(prface.squared_length());
          // Compute afterMach
          double m1 = acos(mu * prface);
          double b1 = obliqueShockAngle(local_mach, m1, gamma);
          double tmach = afterMach(local_mach, b1, m1, gamma);
          // Compute mach wave
          double tmu = asin(1.0 / tmach);
          // Compute plane
          Vector nxz = CGAL::cross_product(mu, prface);
          nxz = nxz / sqrt(nxz.squared_length());
          Vector nc1 = mu;
          Vector nc2 = CGAL::cross_product(nxz, nc1);
          nc2 = nc2 / sqrt(nc2.squared_length());
          Plane fcp(ppair.first, (-sin(tmu+m1) * nc1 + cos(tmu+m1) * nc2));
          // Compute newMach
          Vector nface = streams[i][1]-streams[i][0];
          nface = nface / sqrt(nface.squared_length());
          m1 = acos(mu * nface);
          b1 = obliqueShockAngle(local_mach, m1, gamma);
          tmach = afterMach(local_mach, b1, m1, gamma);
          // Compute mach wave
          tmu = asin(1.0 / tmach);
          // Compute Ray
          nxz = CGAL::cross_product(mu, nface);
          nxz = nxz / sqrt(nxz.squared_length());
          nc1 = mu;
          nc2 = CGAL::cross_product(nxz, nc1);
          nc2 = nc2 / sqrt(nc2.squared_length());
          Ray ff(streams[i][0], (cos(tmu+m1) * nc1 + sin(tmu+m1) * nc2));
          // Compute intersection
          auto apx = intersection(fcp, ff);
          if (apx) {
            start = boost::get<Point>(apx.value());
          } else {
            start = streams[i][0];
          }
        } catch(...) {
          start = streams[i][0];
        }
      }

      double prev_angle = M_PI;
      if (streams[i].size() > 1) {
        //Declare stuff
        double tmu;
        Vector p_n;
        double tba;
        double tmach;
        Vector nxz;
        double tbb;
        Transf rot;
        Vector n;
        Plane p_p;
        double kons = 1;
        unsigned int nst = 0;

        Point apex;
        double ellip;

        //First facet
        Vector tg = streams[i][1] - streams[i][0];
        tg = tg / sqrt(tg.squared_length());
        Vector p_tg = tg;

        double theta = acos(mu * tg);

        if (theta > ma) {
          // std::cout << "Detached shock" << std::endl;

          double maxK;
          double minK;

          std::pair<std::array<double, 2>, std::vector<Vector> >  kpair;
          try {
            kpair = ks.at(streams[i][0]);
          } catch(...) {
            std::array<double, 2> empt = {1.e30, 0};
            std::vector<Vector> tmpv;
            tmpv.push_back(tg);
            tmpv.push_back(CGAL::cross_product(mu, tg));
            tmpv.push_back(mu);

            kpair = std::make_pair(empt, tmpv);
          }
          std::array<double, 2> curvs = kpair.first;
          if (abs(curvs[0]) > abs(curvs[1])) {
            maxK = 1/abs(curvs[1]);
            minK = 1/abs(curvs[0]);
          } else {
            maxK = 1/abs(curvs[0]);
            minK = 1/abs(curvs[1]);
          }

          double sk = minK * 0.143 * exp(3.24 / pow(local_mach,2));
          double ck = minK * 0.386 * exp(4.67 / pow(local_mach,2));

          double rrs = minK * 1.143 * exp(0.54 / pow(local_mach - 1, 1.2));
          double rrc = minK * 1.386 * exp(1.8 / pow(local_mach - 1, 0.75));

          double d = ck * exp(log(sk/ck) / pow(maxK/minK, 0.57987417201950564));
          double rc1 = rrc * exp(log(rrs/rrc) / pow(maxK/minK, 0.57987417201950564));
          double rc2 = 0.9583619656145238 * maxK/minK + rrs - 0.9583619656145238;


          unsigned int j_counter;
          for (j_counter = 1; j_counter < streams[i].size()-1; j_counter++) {
            Vector fv = streams[i][j_counter+1] - streams[i][j_counter];
            fv = fv / sqrt(fv.squared_length());

            double th = acos(mu * fv);

            if (th <= ma) {
              break;
            }
          }

          if (j_counter == streams[i].size()-1) {
            continue;
          }

          nst = j_counter;

          double rc;
          if(model == 1) {
            std::vector<Vector> tried = kpair.second;

            Vector v1 =  streams[i][nst] - streams[i][0];
            v1 = v1 / sqrt(v1.squared_length());

            Vector d1 = (v1 * tried[0]) * tried[0] + (v1 * tried[1]) * tried[1];
            d1 = d1 / sqrt(d1.squared_length());

            double kt = 1/rc2 * pow(d1 * tried[0], 2) +  1/rc1 * pow(d1 * tried[1], 2);
            rc = 1/kt;
          } else if (model == 0) {
            rc = rrs;
          } else {
            rc = rrc;
          }

          Point cs;
          cs = streams[i][0] + (rc - d) * mu;

          // Compute correction factor
          kons = k_factor(local_mach, gamma);

          Point pta = streams[i][nst];
          Point ptb = streams[i][nst + 1];

          Vector fv1 = ptb - pta;
          fv1 = fv1 / sqrt(fv1.squared_length());
          p_tg = fv1;

          double m1 = acos(mu * fv1);

          double b1 = obliqueShockAngle(local_mach, m1, gamma);

          nxz = CGAL::cross_product(mu, pta - streams[i][0]);
          nxz = nxz / sqrt(nxz.squared_length());

          tbb = M_PI / 2 - b1;

          Vector nc1 = mu;
          Vector nc2 = CGAL::cross_product(nxz, nc1);
          nc2 = nc2 / sqrt(nc2.squared_length());

          Point s1 = cs + rc * (-cos(M_PI / 2 - b1) * nc1 + sin(M_PI / 2 - b1) * nc2);

          tmach = afterMach(local_mach, b1, m1, gamma);
          tmu = asin(1.0 / tmach);
          tmu = tmu * kons;

          // Insert first Plane
          n = -sin(b1) * nc1 +  cos(b1) * nc2;

          p_n = nxz;
          p_p = Plane(s1, n);

          double dang = tlen / rc;
          for (double ang = 0; ang < M_PI / 2 - b1; ang += dang) {
            Point tmpt = cs + rc * (-cos(ang) * nc1 + sin(ang) * nc2);
            points.push_back(tmpt+Vector(distribution(generator), distribution(generator), distribution(generator)));
          }
        } else {
          nxz = CGAL::cross_product(mu, tg);
          nxz = nxz / sqrt(nxz.squared_length());
          p_n = nxz;

          double cas = 1;
          if (model == 1) {
            apex = streams[i][0];
            std::pair<std::array<double, 2>, std::vector<Vector> > kval = get_interp(streams[i][1], &tree, kmap);

            Vector nc = CGAL::cross_product(mu, nxz);
            nc = nc / sqrt(nc.squared_length());

            Vector d1 = (nxz * kval.second[0]) * kval.second[0] + (nxz * kval.second[1]) * kval.second[1];
            d1 = d1 / sqrt(d1.squared_length());

            double kn = kval.first[0] * pow(d1 * kval.second[0], 2) +  kval.first[1] * pow(d1 * kval.second[1], 2);
            double K = abs(kn / (nc * kval.second[2]));
            double r2 = abs((apex - streams[i][1]) * nc);

            K = std::min(K,1000.0);

            double a;
            double b;
            if (1/K > r2) {
              a = 1/K;
              b = r2;
              cas = 1;
            } else {
              a = r2;
              b = 1/K;
              cas = 0;
            }

            ellip = sqrt(1 - pow(b,2)/pow(a,2));
          } else if (model == 0){
            ellip = 0;
          } else {
            ellip = 1;
          }

          tba = ellipticShockAngle(local_mach, ellip, theta, gamma, cas);
          prev_angle = tba;
          tmach = ellipticAfterMach(local_mach, ellip, theta, gamma, cas);
          tmu = asin(1.0 / tmach);

          tbb = tba - theta + M_PI/2;

          rot = Transf( cos(tbb) + nxz.x() * nxz.x() * (1 - cos(tbb)), nxz.x() * nxz.y() * (1 - cos(tbb)) - nxz.z() * sin(tbb), nxz.x() * nxz.z() * (1 - cos(tbb)) + nxz.y() * sin(tbb),
                        nxz.y() * nxz.x() * (1 - cos(tbb)) + nxz.z() * sin(tbb), cos(tbb) + nxz.y() * nxz.y() * (1 - cos(tbb)), nxz.y() * nxz.z() * (1 - cos(tbb)) - nxz.x() * sin(tbb),
                        nxz.z() * nxz.x() * (1 - cos(tbb)) - nxz.y() * sin(tbb), nxz.z() * nxz.y() * (1 - cos(tbb)) + nxz.x() * sin(tbb), cos(tbb) + nxz.z() * nxz.z() * (1 - cos(tbb))
                        );

          n = rot(tg);

          p_p = Plane(start, n);
          points.push_back(start);
        }

        for (unsigned int j = nst+1; j < streams[i].size(); j++) {
          ellip = 0;
          int cas = 1;
          if (model == 1) {
            Vector ttg = streams[i][j+1] - streams[i][j];
            ttg = ttg / sqrt(ttg.squared_length());

            Vector nct = CGAL::cross_product(nxz, ttg);
            nct = nct / sqrt(nct.squared_length());
            Plane fcp(streams[i][j+1], nct);

            Ray ff(streams[i][0], -mu);

            auto apx = intersection(fcp, ff);
            if (apx) {
              apex = boost::get<Point>(apx.value());
              std::pair<std::array<double, 2>, std::vector<Vector> > kval = get_interp( streams[i][j+1], &tree, kmap);

              Vector nc = CGAL::cross_product(mu, nxz);
              nc = nc / sqrt(nc.squared_length());

              Vector d1 = (nxz * kval.second[0]) * kval.second[0] + (nxz * kval.second[1]) * kval.second[1];
              d1 = d1 / sqrt(d1.squared_length());

              double kn = kval.first[0] * pow(d1 * kval.second[0], 2) +  kval.first[1] * pow(d1 * kval.second[1], 2);
              double K = abs(kn / (nc * kval.second[2]));
              double r2 = abs((apex - streams[i][j+1]) * nc);

              K = std::min(K,1000.0);

              double a;
              double b;
              if (1/K > r2) {
                a = 1/K;
                b = r2;
                cas = 1;
              } else {
                a = r2;
                b = 1/K;
                cas = 0;
              }

              ellip = sqrt(1 - pow(b,2)/pow(a,2));

            } else {
              std::cout << "Equivalent cone not established on interior point" << std::endl;
              break;
            }
          } else if (model == 0){
            ellip = 0;
          } else {
            ellip = 1;
          }

          Transf rot_1( cos(tmu) + p_n.x() * p_n.x() * (1 - cos(tmu)), p_n.x() * p_n.y() * (1 - cos(tmu)) - p_n.z() * sin(tmu), p_n.x() * p_n.z() * (1 - cos(tmu)) + p_n.y() * sin(tmu),
                        p_n.y() * p_n.x() * (1 - cos(tmu)) + p_n.z() * sin(tmu), cos(tmu) + p_n.y() * p_n.y() * (1 - cos(tmu)), p_n.y() * p_n.z() * (1 - cos(tmu)) - p_n.x() * sin(tmu),
                        p_n.z() * p_n.x() * (1 - cos(tmu)) - p_n.y() * sin(tmu), p_n.z() * p_n.y() * (1 - cos(tmu)) + p_n.x() * sin(tmu), cos(tmu) + p_n.z() * p_n.z() * (1 - cos(tmu))
                        );

          Vector mw = rot_1(p_tg);

          tg = streams[i][j+1] - streams[i][j];
          tg = tg / sqrt(tg.squared_length());
          p_tg = tg;

          theta = acos(mu * tg);

          tba = ellipticShockAngle(local_mach, ellip, theta, gamma, cas);
          if (tba > prev_angle) {
            tba = prev_angle;
          } else {
            prev_angle = tba;
          }
          tmach = ellipticAfterMach(local_mach, ellip, theta, gamma, cas);
          tmu = asin(1.0 / tmach) * kons; // Corrected

          nxz = CGAL::cross_product(tg, mw);
          nxz = nxz / sqrt(nxz.squared_length());
          p_n = nxz;

          tbb = tba - theta + M_PI/2;

          rot = Transf( cos(tbb) + nxz.x() * nxz.x() * (1 - cos(tbb)), nxz.x() * nxz.y() * (1 - cos(tbb)) - nxz.z() * sin(tbb), nxz.x() * nxz.z() * (1 - cos(tbb)) + nxz.y() * sin(tbb),
                        nxz.y() * nxz.x() * (1 - cos(tbb)) + nxz.z() * sin(tbb), cos(tbb) + nxz.y() * nxz.y() * (1 - cos(tbb)), nxz.y() * nxz.z() * (1 - cos(tbb)) - nxz.x() * sin(tbb),
                        nxz.z() * nxz.x() * (1 - cos(tbb)) - nxz.y() * sin(tbb), nxz.z() * nxz.y() * (1 - cos(tbb)) + nxz.x() * sin(tbb), cos(tbb) + nxz.z() * nxz.z() * (1 - cos(tbb))
                        );

          n = rot(tg);

          Ray d(streams[i][j], streams[i][j]+mw);
          auto p1 = intersection(p_p, d);

          if (p1) {
            Point pt = boost::get<Point>(p1.value());
            if ( sqrt((streams[i][0] - pt).squared_length()) < 2*d3) {
              p_p = Plane(pt,n);
              points.push_back(pt+Vector(distribution(generator), distribution(generator), distribution(generator)));
            } else{
              break;
            }
          } else {
            // std::cout << "Streamline following failed" << std::endl;
            break;
          }

        }
      }
    }

    // Filter duplicates
    std::vector<Point>::iterator it;
    it = std::unique (points.begin(), points.end(), pcomp);
    points.resize( std::distance(points.begin(),it) );

    // Mesh
    std::vector<Face> facets;
    double radius_ratio_bound = 2;
    double beta = M_PI/6;
    CGAL::advancing_front_surface_reconstruction(points.begin(), points.end(), std::back_inserter(facets), radius_ratio_bound, beta);

    // Convert to Polyhedron
    Polyhedron St;
    BuildPoly<HalfedgeDS> builder(&points, &facets);
    St.delegate( builder);
    St.keep_largest_connected_components(1);
    S.push_back(St);

    Polyhedron Sst;
    CGAL::convex_hull_3(points.begin(), points.end(), Sst);
    Sch.push_back(Sst);

    // Mean after shock mach
    double tamach = 0;
    int nm = 0;
    for ( Facet_iterator f = St.facets_begin(); f != St.facets_end(); ++f ) {
      Halfedge_handle h = f->halfedge();

      Vector norm;
      try {
        norm = CGAL::unit_normal(h->vertex()->point(),
                                 h->next()->vertex()->point(),
                                 h->prev()->vertex()->point());
      } catch(...) {
        continue;
      }

      double b = M_PI/2 - acos(-mu * norm);

      double t =  atan((2/tan(b) * (pow(local_mach, 2) * pow(sin(b), 2)-1)) / (pow(local_mach, 2) * (gamma + cos(2*b))+2));
      double amach = 1.0 / sin(b - t) * sqrt((1.0 + (gamma - 1.0) / 2.0 * pow(local_mach, 2) * pow(sin(b), 2)) / (gamma * pow(local_mach, 2) * pow(sin(b), 2) - (gamma - 1.0) / 2.0));
      amach = abs(amach);
      if (amach > 1) {
        tamach += amach; // Change to area averaged
        ++nm;
      }
    }

    after_mach.push_back(tamach/nm);
  }

  // Write output
  return S;
}
