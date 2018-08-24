#include "tmesh.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <fstream>



namespace T_MESH {


#define TVI1(a) (TMESH_TO_INT(((Triangle *)a->data)->v1()->x))
#define TVI2(a) (TMESH_TO_INT(((Triangle *)a->data)->v2()->x))
#define TVI3(a) (TMESH_TO_INT(((Triangle *)a->data)->v3()->x))


    void setVertexDouble(double *vertices, float x, float y, float z, int num_vertex, int idx) {
        vertices[idx + 0*num_vertex] = (double) x;
        vertices[idx + 1*num_vertex] = (double) y;
        vertices[idx + 2*num_vertex] = (double) z;
//        printf("v1:%f\t v2:%f\t v3:%f\n", x, y, z);
    }

    void setFacesDouble(double *faces, int x, int y, int z, int num_faces, int idx) {
        faces[idx + 0*num_faces] = (double) x;
        faces[idx + 1*num_faces] = (double) y;
        faces[idx + 2*num_faces] = (double) z;
//        printf("f1:%d\t f2:%d\t f3:%d\n", x, y, z);
    }

    int Basic_TMesh::loadDouble(
            const double *in_vertices,
            const size_t n_vertex,
            const double *in_faces,
            const size_t n_faces,
            const bool doupdate) {

        Node *n;
        float x, y, z;
        size_t i = 0;
        int i1, i2, i3, triangulate = 0;
        Vertex *v;

        if (n_vertex < 3) TMesh::error("\nloadEigen: Can't load objects with less than 3 vertices.\n");
        if (n_faces < 1) TMesh::error("\nloadEigen: Can't load objects with no faces.\n");

        for (i = 0; i < n_vertex; i++) {
            x = (float) in_vertices[i];
            y = (float) in_vertices[i + n_vertex];
            z = (float) in_vertices[i + 2 * n_vertex];
            V.appendTail(newVertex(x, y, z));
        }

        ExtVertex **var = (ExtVertex **) malloc(sizeof(ExtVertex * ) * n_vertex);
        i = 0; FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);

        TMesh::begin_progress();
        for (i = 0; i < n_faces; i++) {
            i1 = (int) in_faces[i];
            i2 = (int) in_faces[i + n_faces];
            i3 = (int) in_faces[i + 2 * n_faces];

            if ((i % 1000) == 0) TMesh::report_progress("Loading ..%d%%", (i * 100) / (n_vertex * 2));

            if (i1 < 0 || i2 < 0 || i3 < 0 || i1 > (n_vertex - 1) || i2 > (n_vertex - 1) || i3 > (n_vertex - 1))
                TMesh::error("\nloadDouble: Invalid index at face %d!\n i1:%d\t i2:%d\t i3:%d\t nV:%d\n", i, i1, i2, i3, n_vertex);
            if (i1 == i2 || i2 == i3 || i3 == i1)
                TMesh::warning("\nloadDouble: Coincident indexes at triangle %d! Skipping.\n", i);
            else if (!CreateIndexedTriangle(var, i1, i2, i3))
                TMesh::warning("\nloadDouble: This shouldn't happen!!! Skipping triangle.\n");

        }

        TMesh::end_progress();

        closeLoading(i, var, (triangulate != 0));

        if (doupdate) eulerUpdate();

        return 0;
    }


    void Basic_TMesh::closeLoading(int loaded_faces, ExtVertex **var, bool triangulate) {
        int i, n_vertex = V.numels();

        if (var != NULL) {
            for (i = 0; i < n_vertex; i++) delete (var[i]);
            free(var);
        }

        if (loaded_faces) {
            TMesh::info("Loaded %d vertices and %d faces.\n", n_vertex, loaded_faces);
            if (triangulate) TMesh::warning("Some polygonal faces needed to be triangulated.\n");
            fixConnectivity();
        }

        d_boundaries = d_handles = d_shells = 1;
    }

//    <int Basic_TMesh::getNumVertices() {
//        return V.numels();
//    }
//
//    int Basic_TMesh::getNumFaces() {
//        return T.numels();
//    }


    int Basic_TMesh::exportDouble(
            double *out_vertices,
            double *out_faces) {
        int i;
        int idx;
        Node *n;
        coord *ocds;
        Vertex *v;

        int num_vertex = V.numels();
        int num_faces = T.numels();

        idx = 0;
        FOREACHVERTEX(v, n)
        setVertexDouble(out_vertices,TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z), num_vertex, idx++);

        ocds = new coord[V.numels()];
        i = 0;
        FOREACHVERTEX(v, n)
        ocds[i++] = v->x;
        i = 0;
        FOREACHVERTEX(v, n)
        v->x = i++;

        idx = 0;
        FOREACHNODE(T, n)
        setFacesDouble(out_faces, TVI1(n), TVI2(n), TVI3(n), num_faces, idx++);

        i = 0;
        FOREACHVERTEX(v, n)
        v->x = ocds[i++];
        delete[] ocds;

        return 0;
    }

} //namespace T_MESH
