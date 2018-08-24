//
// Created by david on 13.08.18.
//

#include "tmesh.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <mex.h>

using namespace T_MESH;


double closestPair(List *bl1, List *bl2, Vertex **closest_on_bl1, Vertex **closest_on_bl2) {
    Node *n, *m;
    Vertex *v, *w;
    double adist, mindist = DBL_MAX;

    FOREACHVVVERTEX(bl1, v, n)FOREACHVVVERTEX(bl2, w, m)if ((adist = w->squaredDistance(v)) < mindist) {
                mindist = adist;
                *closest_on_bl1 = v;
                *closest_on_bl2 = w;
            }

    return mindist;
}

bool joinClosestComponents(Basic_TMesh *tin) {
    Vertex *v, *w, *gv, *gw;
    Triangle *t, *s;
    Node *n;
    List triList, boundary_loops, *one_loop;
    List **bloops_array;
    int i, j, numloops;

    // Mark triangles with connected component's unique ID
    i = 0;
    FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
    FOREACHVTTRIANGLE((&(tin->T)), t, n) if (t->info == NULL) {
            i++;
            triList.appendHead(t);
            t->info = (void *) i;

            while (triList.numels()) {
                t = (Triangle *) triList.popHead();
                if ((s = t->t1()) != NULL && s->info == NULL) {
                    triList.appendHead(s);
                    s->info = (void *) i;
                }
                if ((s = t->t2()) != NULL && s->info == NULL) {
                    triList.appendHead(s);
                    s->info = (void *) i;
                }
                if ((s = t->t3()) != NULL && s->info == NULL) {
                    triList.appendHead(s);
                    s->info = (void *) i;
                }
            }
        }

    if (i < 2) {
        FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
        //   JMesh::info("Mesh is a single component. Nothing done.");
        return false;
    }

    FOREACHVTTRIANGLE((&(tin->T)), t, n) {
        t->v1()->info = t->v2()->info = t->v3()->info = t->info;
    }

    FOREACHVVVERTEX((&(tin->V)), v, n) if (!IS_VISITED2(v) && v->isOnBoundary()) {
            w = v;
            one_loop = new List;
            do {
                one_loop->appendHead(w);
                MARK_VISIT2(w);
                w = w->nextOnBoundary();
            } while (w != v);
            boundary_loops.appendHead(one_loop);
        }
    FOREACHVVVERTEX((&(tin->V)), v, n) UNMARK_VISIT2(v);

    bloops_array = (List **) boundary_loops.toArray();
    numloops = boundary_loops.numels();

    int numtris = tin->T.numels();
    double adist, mindist = DBL_MAX;

    gv = NULL;
    for (i = 0; i < numloops; i++)
        for (j = 0; j < numloops; j++)
            if (((Vertex *) bloops_array[i]->head()->data)->info != ((Vertex *) bloops_array[j]->head()->data)->info) {
                adist = closestPair(bloops_array[i], bloops_array[j], &v, &w);
                if (adist < mindist) {
                    mindist = adist;
                    gv = v;
                    gw = w;
                }
            }

    if (gv != NULL) tin->joinBoundaryLoops(gv, gw, 1, 0);

    FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
    FOREACHVVVERTEX((&(tin->V)), v, n) v->info = NULL;

    free(bloops_array);
    while ((one_loop = (List *) boundary_loops.popHead()) != NULL) delete one_loop;

    return (gv != NULL);
}


void mexFunction(
        int nlhs,
        mxArray *plhs[],
        int nrhs,
        const mxArray *prhs[]
) {
    double *V_in;
    double *F_in;
    double *V_out;
    double *F_out;
    size_t num_vertices_in, num_faces_in;
    size_t num_vertices_out, num_faces_out;
//    bool join_multiple_components = false;
    bool auto_update = true;

    if (mxGetN(prhs[0]) != 3) mexErrMsgTxt("V_in must be #V_in by 3");
    if (mxGetN(prhs[1]) != 3) mexErrMsgTxt("F_in must be #F_in by 3");
    V_in = mxGetDoubles(prhs[0]);
    F_in = mxGetDoubles(prhs[1]);
    num_vertices_in = mxGetM(prhs[0]);
    num_faces_in = mxGetM(prhs[1]);
    mexPrintf("Parsing inputs completed\n");

    TMesh::init(); // This is mandatory
    TMesh::app_name = "mat_meshfix";
    TMesh::app_version = "1.0";
    TMesh::app_year = "2018";
    TMesh::app_authors = "David Dreher";
    TMesh::app_maillist = "david.dreher@rocketmail.com";

    mexPrintf("TMesh init completed\n");

    clock_t beginning = clock();

    Basic_TMesh tin;
    mexPrintf("TMesh object constructed\n");

    if (tin.loadDouble(V_in, num_vertices_in, F_in, num_faces_in, auto_update) != 0)
        TMesh::error("Can't import mesh.\n");
    mexPrintf("Loaded double arrays\n");

//    if (join_multiple_components) {
//        TMesh::info("\nJoining input components ...\n");
//        TMesh::begin_progress();
//        while (joinClosestComponents(&tin)) TMesh::report_progress("Num. components: %d       ", tin.shells());
//        TMesh::end_progress();
//        tin.deselectTriangles();
//    }

    // Keep only the largest component (i.e. with most triangles)
    int sc = tin.removeSmallestComponents();
    if (sc) TMesh::warning("Removed %d small components\n", sc);

    // Fill holes
    if (tin.boundaries()) {
        TMesh::warning("Patching holes\n");
        tin.fillSmallBoundaries(0, true);
    }

    // Run geometry correction
    if (!tin.boundaries()) TMesh::warning("Fixing degeneracies and intersections...\n");
    if (tin.boundaries() || !tin.meshclean()) TMesh::warning("MeshFix could not fix everything.\n", sc);

    mexPrintf("Number vertices: %d\n", tin.V.numels() );
    mexPrintf("Number faces: %d\n", tin.T.numels() );

    num_vertices_out = (size_t) tin.V.numels();
    num_faces_out = (size_t) tin.T.numels();
    plhs[0] = mxCreateDoubleMatrix(num_vertices_out, 3, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(num_faces_out, 3, mxREAL);

    V_out = mxGetDoubles(plhs[0]);
    F_out = mxGetDoubles(plhs[1]);
    mexPrintf("Exporting Manifold Oriented Triangulation ...\n");
    if (tin.exportDouble(V_out, F_out) != 0)
        TMesh::error("Can't export mesh.\n");

//
//
//    mxSetDoubles(plhs[0], V_out);
//    mxSetDoubles(plhs[1], F_out);

    printf("Elapsed time: %d ms\n", clock() - beginning);

//    int i;
//    mexPrintf("New vertices(%d):\n", num_vertices_out);
//    for (i = 0; i < num_vertices_out; i++) {
//        mexPrintf("%f\t%f\t%f\t\n", V_out[i], V_out[i+num_vertices_out], V_out[i+2*num_vertices_out]);
//    }
//
//    mexPrintf("New faces(%d):\n", num_faces_out);
//    for (i = 0; i < num_faces_out; i++) {
//        mexPrintf("%d\t%d\t%d\t\n", (int) F_out[i], (int) F_out[i+num_faces_out], (int) F_out[i+2*num_faces_out]);
//    }

}
