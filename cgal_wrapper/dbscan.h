// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include "./types.h"
#include "kdtree.h"

void expandCluster(unsigned int pt_index, std::vector<std::pair<size_t,double>> & NeighborPts, int C, double eps, unsigned int MinPts, std::vector<std::pair<Point, int>> & D, KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud>,PointCloud,3> & treeSI) {
  D[pt_index].second = C;
  while (!NeighborPts.empty()) {
    std::pair<size_t,double> tmpp = NeighborPts.back();
    NeighborPts.pop_back();
    size_t tindx = tmpp.first;

    if(D[tindx].second ==-2) {
      D[tindx].second = C;
      Point qp = D[tindx].first;
      std::vector<double> query_pt = {qp.x(), qp.y(), qp.z()};
      std::vector<std::pair<size_t,double>>   NeighborPts2;
      nanoflann::SearchParams params;

      size_t nMatches = treeSI.radiusSearch(&query_pt[0], eps*eps, NeighborPts2, params);
      nMatches += 1;
      NeighborPts2.push_back(std::pair<size_t, double>(tindx, 0.0));

      if(nMatches >= MinPts) {
        NeighborPts.insert(NeighborPts.end(), NeighborPts2.begin(), NeighborPts2.end());
      }
    }
    if(D[tindx].second ==-1) {
      D[tindx].second = C;
    }
  }
}

int DBSCAN(std::vector<std::pair<Point, int>> &D, double eps, unsigned int MinPts) {
  PointCloud cloud;
  cloud.pts.resize(D.size());

  for (unsigned int i = 0; i < D.size(); i++) {
    cloud.pts[i] = D[i].first;
  }

  KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud>,PointCloud,3> treeSI(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  treeSI.buildIndex();

  nanoflann::SearchParams params;

  int C = 0;
  for(unsigned int i=0; i< D.size(); i++) {
    if (D[i].second != -2) {
      continue;
    }

    Point qp = D[i].first;
    std::vector<double> query_pt = {qp.x(), qp.y(), qp.z()};
    std::vector<std::pair<size_t,double>>   NeighborPts;

    size_t nMatches = treeSI.radiusSearch(&query_pt[0], eps*eps, NeighborPts, params);
    nMatches += 1;
    NeighborPts.push_back(std::pair<size_t, double>(i, 0.0));

    if(nMatches < MinPts) {
      D[i].second = -1;
    }
    else {
      C += 1;
      expandCluster( i, NeighborPts, C, eps, MinPts, D, treeSI);
    }
  }

  return C;
}

#endif /* DBSCAN_H */
