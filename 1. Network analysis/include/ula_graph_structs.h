/************************************************************
 *			Struct definition
 * Note: It is compatible with Pol Colomer's RandNetGen definition for compatibility purposes
 * To make them fully compatible just take out commented code
 *************************************************************/


#ifndef ULA_GRAPH_STRUCTS_H
#define ULA_GRAPH_STRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//----- A NODE --------//
typedef struct NODE{
	
	int idnum;
	int kin;
	int kout;
	int sin;
	int sout;
	int *out;
	int *in;
	int mem_in;
	int mem_out;
	int* w_out;
	int* w_in;
	//// 	For spatial networks if needed
	//double x;
	//double y;
	//double loc_x;
	//double loc_y;
	//// 	Compatibility with NetRandGen
	//int tri;    /// number of triangles or common neighbours
	

}NODE;

//----- AN EDGE --------//
typedef struct EDGE{
	
	int s;	/// source
	int d;	/// destination
	int w;	/// weight
}EDGE;

//----- A GRAPH --------//
typedef struct W_GRAPH{
	
	int N_nodes;				/// number of nodes
	NODE* node;			/// vector with the nodes
    int E;				/// number of edges
	EDGE* edge;			/// vector with all the edges
	//// 	Compatibility with NetRandGen
	//int linksWrepeat;	/// number of links with repetitions
	//int max_k;			/// maximum degree
    //int loops;			/// number of self edges
	//int* pk;            /// degree distribution
    //double* ck;         /// clustering spectrum
    //double  Ccoef;      /// clustering coefficient
    //double triangles;   /// number of triangles divided by N
    //double* Knn;        /// the average neighbour degree

}W_GRAPH;

#endif
