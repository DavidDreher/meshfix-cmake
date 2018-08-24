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


    void setVertexDouble(double *vertices, float x, float y, float z, int num_vertex) {
        vertices[0] = (double) x;
        vertices[num_vertex] = (double) y;
        vertices[2 * num_vertex] = (double) z;
    }

    void setFacesDouble(double *faces, int x, int y, int z, int num_faces) {
        faces[0] = (double) x;
        faces[num_faces] = (double) y;
        faces[2 * num_faces] = (double) z;
    }

    int Basic_TMesh::loadDouble(
            const double *in_vertices,
            const int n_vertex,
            const double *in_faces,
            const int n_faces,
            const bool doupdate) {

        Node *n;
        float x, y, z;
        int i, i1, i2, i3, triangulate = 0;
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
                TMesh::error("\nloadDouble: Invalid index at face %d!\n", i);
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


    int Basic_TMesh::exportDouble(
            double *out_vertices,
            int *num_vertex,
            double *out_faces,
            int *num_faces) {
        int i;
        Node *n;
        coord *ocds;
        Vertex *v;

        *num_vertex = (int) V.numels();
        *num_faces = (int) T.numels();

        out_vertices = (double *) malloc(sizeof(double) * *num_vertex * 3);
        FOREACHVERTEX(v, n)
        setVertexDouble(out_vertices,TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z), *num_vertex);

        ocds = new coord[V.numels()];
        i = 0;
        FOREACHVERTEX(v, n)
        ocds[i++] = v->x;
        i = 0;
        FOREACHVERTEX(v, n)
        v->x = i++;

        out_faces = (double *) malloc(sizeof(double) * *num_faces * 3);

        FOREACHNODE(T, n)
        setFacesDouble(out_faces, TVI1(n), TVI2(n), TVI3(n), *num_faces);

        i = 0;
        FOREACHVERTEX(v, n)
        v->x = ocds[i++];
        delete[] ocds;

        return 0;
    }

} //namespace T_MESH
