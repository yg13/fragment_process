// This is brcv/shp/dbsk2d/pro/dbsk2d_compute_containment_graph_process.h
#ifndef dbsk2d_compute_containment_graph_process_h_
#define dbsk2d_compute_containment_graph_process_h_

//:
// \file
// \brief This process computes a containment graph data structure from
//        a vsol2d storage
//
// \author Maruthi Narayanan
// \date 09/21/12
//
// \verbatim
//  Modifications
//
// \endverbatim

#include "../../bpro1/bpro1_process.h"
#include "../../bpro1/bpro1_parameters.h"
#include "../dbsk2d_ishock_graph_sptr.h"
#include "../dbsk2d_bnd_contour_sptr.h"
#include "../pro/dbsk2d_shock_storage_sptr.h"
#include <vgl/vgl_polygon.h>

class dbsk2d_ishock_bpoint;

class dbsk2d_compute_containment_graph_process : public bpro1_process 
{

public:
  //: Constructor
  dbsk2d_compute_containment_graph_process();
  
  //: Destructor
  virtual ~dbsk2d_compute_containment_graph_process();

  //: Clone the process
  virtual bpro1_process* clone() const;

  vcl_string name();

  vcl_vector< vcl_string > get_input_type();
  vcl_vector< vcl_string > get_output_type();

  int input_frames();
  int output_frames();

  bool execute();
  bool finish();
  
private:

  void debug_cost(dbsk2d_ishock_graph_sptr ishock_graph,
                  vcl_string filename);

  void debug_frags(dbsk2d_ishock_graph_sptr ishock_graph,
                   vcl_string filename);

  void pre_process_contours(dbsk2d_ishock_graph_sptr ishock_graph,
                            double preprocess_threshold,
                            double gap_distance,
                            bool remove_closed,
                            bool train);
  
  void create_dependency_graph(dbsk2d_ishock_graph_sptr ishock_graph,
                               vcl_string filename);

  void pre_process_gap4(dbsk2d_shock_storage_sptr output_shock);

  dbsk2d_bnd_contour_sptr get_contour(dbsk2d_ishock_bpoint* bp);

  void region_stats(vgl_polygon<double>& poly,vcl_vector<double>& stats,
                    dbsk2d_ishock_graph_sptr ishock_graph);

};

#endif  //dbsk2d_compute_containment_graph_process_h_
