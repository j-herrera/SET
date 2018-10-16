import numpy as np

def read_ply(in_ply):
    s_file = open(in_ply, "r")

    line = s_file.readline()
    if (line.strip() != "ply"):
        raise RuntimeError("Not a ply file")

    line = s_file.readline()
    if (line.strip() != "format ascii 1.0"):
        raise RuntimeError("Only ascii accepted")

    while "element vertex" not in line:
        line = s_file.readline()

    nv = int(line[14:])

    line = s_file.readline()
    while "element face" not in line:
        line = s_file.readline()

    nf = int(line[12:])

    line = s_file.readline()
    while (line.strip() != "end_header"):
        line = s_file.readline()

    verts = np.zeros([nv,3])
    # assumes properties are x y z ...
    for i in range(nv):
        line = s_file.readline()
        verts[i,:] = np.fromstring(line, count = 3, sep=' ')

    faces = np.zeros([nf,3], dtype='int')
    for i in range(nf):
        line = s_file.readline()
        data = np.fromstring(line, dtype = 'int', sep=' ')
        if data[0] != 3:
            raise RuntimeError("Only triangular meshes accepted")
        else:
            faces[i,:] = data[1:]

    s_file.close()

    return {'vertices': verts, 'faces':faces}

def write_ply(out_ply, data):
    s_file = open(out_ply, "w")

    s_file.write("ply\n")
    s_file.write("format ascii 1.0\n")

    s_file.write("element vertex " + str(data[1].shape[0]) + "\n")
    s_file.write("property float x\n")
    s_file.write("property float y\n")
    s_file.write("property float z\n")

    s_file.write("element face " + str(data[0].shape[0]) + "\n")
    s_file.write("property list uchar uint vertex_indices\n")
    s_file.write("end_header\n")

    for i in range(data[1].shape[0]):
        s_file.write(str(data[1][i,:])[1:-1]+"\n")

    for i in range(data[0].shape[0]):
        s_file.write(str(data[0][i,:].size) + " " + np.array2string(data[0][i,:])[1:-1]+"\n")

    s_file.close()
