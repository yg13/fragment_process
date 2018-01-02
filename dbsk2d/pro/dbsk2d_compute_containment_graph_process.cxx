// This is brcv/shp/dbsk2d/pro/dbsk2d_compute_containment_graph_region_process.cxx

//:
// \file

#include "dbsk2d_compute_containment_graph_process.h"
#include "../pro/dbsk2d_shock_storage.h"
#include "../pro/dbsk2d_compute_ishock_process.h"
#include "../dbsk2d_containment_graph.h"
#include "../dbsk2d_transform_manager.h"
#include "../algo/dbsk2d_compute_shocks.h"
#include "../algo/dbsk2d_ishock_gap_detector.h"
#include "../algo/dbsk2d_ishock_transform_sptr.h"
#include "../algo/dbsk2d_ishock_gap_transform.h"
#include "../algo/dbsk2d_ishock_loop_transform.h"
#include "../algo/dbsk2d_ishock_gap4_transform.h"
#include "../algo/dbsk2d_ishock_detector.h"
#include "../dbsk2d_shock_grouped_ishock_edge.h"
#include "../algo/dbsk2d_ishock_grouping_transform.h"
#include "../algo/dbsk2d_bnd_preprocess.h"
#include "../dbsk2d_bnd_utils.h"
#include "../algo/dbsk2d_prune_ishock.h"
#include "../../dbsol/dbsol_interp_curve_2d.h"
#include "../../dbsol/algo/dbsol_curve_algs.h"
#include <vgl/vgl_area.h>
#include <vgl/vgl_line_segment_2d.h>

#include "../../vidpro1/storage/vidpro1_vsol2D_storage_sptr.h"
#include "../../vidpro1/storage/vidpro1_vsol2D_storage.h"
#include "../../vidpro1/storage/vidpro1_image_storage.h"
#include "../../vidpro1/storage/vidpro1_image_storage_sptr.h"
#include <vil/vil_image_resource.h>


#include <vnl/vnl_random.h>
#include <vgl/vgl_lineseg_test.txx>
#include <vgl/vgl_closest_point.h>
#include <vgl/vgl_box_2d.h>
#include <vsol/vsol_point_2d.h>
#include <vsol/vsol_line_2d.h>
#include <vsol/vsol_polyline_2d.h>
#include <vsol/vsol_box_2d.h>
#include <vsol/vsol_polygon_2d.h>
#include <bsol/bsol_algs.h>

#include <vgl/vgl_distance.h>
#include <vgl/algo/vgl_fit_lines_2d.h>
#include <vul/vul_timer.h>

#include "../../dbsol/dbsol_file_io.h"
//: Constructor
dbsk2d_compute_containment_graph_process::
dbsk2d_compute_containment_graph_process()
{
    if (
            !parameters()->add( "Region id:" , "-region_id" , (int) 265 ) ||
            !parameters()->add( "Expand Single:" ,
                                "-expand_single" , (bool) false ) ||
            !parameters()->add( "Add bbox:" ,
                                "-add_bbox", (bool) true) ||
            !parameters()->add( "Expansion Type: Implicit(1) Explicit(0)" ,
                                "-expand_type" , (int) 1 ) ||
            !parameters()->add( "Threshold on Probability of Transform:" ,
                                "-transform_threshold" ,
                                (double) 0.075 ) ||
            !parameters()->add( "Threshold on Probability of Path:" ,
                                "-path_threshold" ,
                                (double) 0.05 ) ||
            !parameters()->add( "Threshold on Preprocessing transforms:" ,
                                "-preprocess_threshold" ,
                                (double) 0.1 ) ||
            !parameters()->add( "Minimum Gap distance:" ,
                                "-gap_distance" ,
                                (double) 3.5 ) ||
            !parameters()->add( "Euler Sprial Completion" , "-ess" ,
                                (double)0.30) ||
            !parameters()->add( "Remove closed contours" , "-closed" ,
                                (bool)false) ||
            !parameters()->add( "Quad" , "-quad" ,
                                (int) -1) ||
            !parameters()->add( "Loop cost type (0,1)" , "-loop_cost" ,
                                (unsigned int) 1) ||
            !parameters()->add( "Expand Outside" , "-outside" ,
                                (bool) false) ||
            !parameters()->add( "Train" , "-train" ,
                                (bool) false) ||
            !parameters()->add( "Debug" , "-debug" ,
                                (bool) false) ||
            !parameters()->add( "Show shock" , "-show_shock" ,
                                (bool) false) ||
            !parameters()->add( "Output folder:" ,
                                "-output_folder", bpro1_filepath("", "")) ||
            !parameters()->add( "Output file prefix:" ,
                                "-output_prefix", vcl_string("")) ||
            !parameters()->add("If grow bbox based on contours", "-grow_by_bnd", (bool)false)

            )
    {
        vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;

    }

}

//: Destructor
dbsk2d_compute_containment_graph_process::
~dbsk2d_compute_containment_graph_process()
{
}

//: Clone the process
bpro1_process*
dbsk2d_compute_containment_graph_process::clone() const
{
    return new dbsk2d_compute_containment_graph_process(*this);
}

vcl_string
dbsk2d_compute_containment_graph_process::name()
{
    return "Compute Containment Graph Process";
}

vcl_vector< vcl_string >
dbsk2d_compute_containment_graph_process::get_input_type()
{
    vcl_vector< vcl_string > to_return;
    to_return.push_back( "vsol2D" );
    to_return.push_back( "image" );
    return to_return;
}

vcl_vector< vcl_string >
dbsk2d_compute_containment_graph_process::get_output_type()
{
    vcl_vector< vcl_string > to_return;
    return to_return;
}

int dbsk2d_compute_containment_graph_process::input_frames()
{
    return 1;
}

int dbsk2d_compute_containment_graph_process::output_frames()
{
    return 1;
}

bool dbsk2d_compute_containment_graph_process::execute()
{

    // Start timer
    vul_timer t;
	clock_t timer;
	timer = clock();

    double transform_threshold(0.5);
    double path_threshold(0.2);
    vcl_string output_prefix;
    bool remove_closed(false);
    bool add_bbox(true);
    double preprocess_threshold(0.12);
    unsigned int loop_cost(1);
    double gap_distance(2.0);
    bool outside(false);
    bool train(false);
    int quad(-1);
    bool debug(false);
    bool show_shock(false);
	bool grow_by_bnd(false);


    bpro1_filepath output_folder_filepath;
    this->parameters()->get_value("-output_folder", output_folder_filepath);
    vcl_string output_folder = output_folder_filepath.path;


    parameters()->get_value( "-transform_threshold" , transform_threshold );
    parameters()->get_value( "-path_threshold" , path_threshold );
    parameters()->get_value( "-output_prefix", output_prefix);
    parameters()->get_value( "-closed" ,         remove_closed );
    parameters()->get_value( "-add_bbox", add_bbox);
    parameters()->get_value( "-preprocess_threshold" , preprocess_threshold );
    parameters()->get_value( "-loop_cost" , loop_cost);
    parameters()->get_value( "-gap_distance" , gap_distance );
    parameters()->get_value( "-outside" ,  outside );
    parameters()->get_value( "-train" , train );
    parameters()->get_value( "-debug" , debug );
    parameters()->get_value( "-show_shock" , show_shock );
    parameters()->get_value( "-quad" , quad );
    parameters()->get_value("-grow_by_bnd", grow_by_bnd);

    bool status = true;

    // 1) get input storage class
    vidpro1_vsol2D_storage_sptr input_vsol;
    input_vsol.vertical_cast(input_data_[0][0]);

    //1) get input storage classes
    vidpro1_image_storage_sptr frame_image;
    frame_image.vertical_cast(input_data_[0][1]);

    if ( !dbsk2d_transform_manager::Instance().gPb_loaded() )
    {
        vcl_cout<<"gPb data not loaded, setting all thresholds to negative"
                <<vcl_endl;

        transform_threshold  = -1.0;
        preprocess_threshold = -1.0;
        path_threshold       = -1.0;
    }

    // Set other thresholds
    dbsk2d_transform_manager::Instance().set_threshold(transform_threshold);
    dbsk2d_transform_manager::Instance().set_output_frag_folder(output_folder);
    dbsk2d_transform_manager::Instance().set_output_prefix(output_prefix);

    vil_image_resource_sptr image(0);
    if ( frame_image)
    {
        image=frame_image->get_image();

        vcl_string output_binary_file = output_folder+"/"+
                                        output_prefix+"_regions_binary.bin";
        vcl_string output_binary_regions_contours_file = output_folder+"/"+
                                                         output_prefix+"_regions_contours_binary.bin";
        vcl_string output_region_stats_file = output_folder+"/"+
                                              output_prefix+"_regions_stats.bin";

        dbsk2d_transform_manager::Instance().set_image(image);
        dbsk2d_transform_manager::Instance().start_binary_file(
                output_binary_file);
        dbsk2d_transform_manager::Instance().start_region_file(
                output_binary_regions_contours_file);
        dbsk2d_transform_manager::Instance().start_region_stats_file(
                output_region_stats_file);

    }

    // Redo input vsol storage

    // 2) Set id equivalent to what they are
    vcl_vector< vsol_spatial_object_2d_sptr > vsol_list =
            input_vsol->all_data();

    // 3) Keep track of objects to remove
    vcl_vector< vsol_spatial_object_2d_sptr > objects_to_erase;
    vcl_map< double, vsol_spatial_object_2d_sptr > low_prob_curves;

    vsol_box_2d_sptr bbox=0;

    // 4) Adding bounding box if needed
    if ( add_bbox )
    {

        // create new bounding box
        bbox = new vsol_box_2d();

        // Grab image to see boundaries
        vil_image_resource_sptr image = frame_image->get_image();

        // Enlarge bounding box from size
        // Calculate xcenter, ycenter
        double xcenter = image->ni()/2.0;
        double ycenter = image->nj()/2.0;

        // Translate to center and scale
        double xmin_scaled = ((0-xcenter)*2)+xcenter;
        double ymin_scaled = ((0-ycenter)*2)+ycenter;
        double xmax_scaled = ((image->ni()-xcenter)*2)+xcenter;
        double ymax_scaled = ((image->nj()-ycenter)*2)+ycenter;

        vcl_cout<<xmin_scaled<<" "<<ymin_scaled<<vcl_endl;
        bbox->add_point(xmin_scaled,ymin_scaled);
        bbox->add_point(xmax_scaled,ymax_scaled);

        vcl_cout << "bbox (minx, miny) (maxx, maxy) (width, height): "
                 << "("   << bbox->get_min_x() << ", " << bbox->get_min_y()
                 << ") (" << bbox->get_max_x() << ", " << bbox->get_max_y()
                 << ") ("
                 << bbox->width() << ", "
                 << bbox->height() << ")"<<vcl_endl;

        // Add to vidpro storage this new bounding box
        vsol_polygon_2d_sptr box_poly = bsol_algs::poly_from_box(bbox);
        input_vsol->add_object(box_poly->cast_to_spatial_object());
	    box_poly->print(vcl_cout);
        bbox=0;
    }

	else if( grow_by_bnd )
	{
        // create new bounding box
        bbox = new vsol_box_2d();

		for (unsigned int b = 0; b < vsol_list.size(); b++) {
			bbox->grow_minmax_bounds(vsol_list[b]->get_bounding_box());
		}

        vil_image_resource_sptr image = frame_image->get_image();
		// Enlarge bounding box from size
		// Calculate xcenter, ycenter
		double xcenter = bbox->width() / 2.0;
		double ycenter = bbox->height() / 2.0;

		// Translate to center and scale
		double xmin_scaled = ((bbox->get_min_x() - xcenter) * 2) + xcenter;
		double ymin_scaled = ((bbox->get_min_y() - ycenter) * 2) + ycenter;
		double xmax_scaled = ((bbox->get_max_x() - xcenter) * 2) + xcenter;
		double ymax_scaled = ((bbox->get_max_y() - ycenter) * 2) + ycenter;

		bbox->add_point(xmin_scaled, ymin_scaled);
		bbox->add_point(xmax_scaled, ymax_scaled);

		vcl_cout << "bbox (minx, miny) (maxx, maxy) (width, height): "
		         << "(" << bbox->get_min_x() << ", " << bbox->get_min_y()
		         << ") (" << bbox->get_max_x() << ", " << bbox->get_max_y()
		         << ") ("
		         << bbox->width() << ", "
		         << bbox->height() << ")" << vcl_endl;

		// Add to vidpro storage this new bounding box
		vsol_polygon_2d_sptr box_poly = bsol_algs::poly_from_box(bbox);
		input_vsol->add_object(box_poly->cast_to_spatial_object());
		box_poly->print(vcl_cout);
        bbox=0;
	}

    double max_gPb_value=0.0;
    vcl_vector<vgl_polygon<double> > closed_polys;

    for ( unsigned int i=0; i < vsol_list.size() ; ++i)
    {
        if( vsol_list[i]->cast_to_curve())
        {
            if( vsol_list[i]->cast_to_curve()->cast_to_polyline() )
            {
                vsol_polyline_2d_sptr curve =
                        vsol_list[i]->cast_to_curve()->cast_to_polyline();

                double prob = dbsk2d_transform_manager::Instance().
                        transform_probability(curve);
                low_prob_curves[prob]=vsol_list[i];

                if ( curve->p0()->get_p() == curve->p1()->get_p())
                {
                    // This represents a closed curve
                    // Lets write it out and not use it

                    // Do not include this in the contour set
                    objects_to_erase.push_back(vsol_list[i]);

                    // Grab all points from vsol
                    vcl_vector<vgl_point_2d<double> > points;
                    for (unsigned p=0; p<curve->size(); p++)
                    {
                        points.push_back( curve->vertex(p)->get_p());
                    }

                    // Create polygon
                    vgl_polygon<double> poly(points);


                    int min_fit_length = 2;
                    vgl_fit_lines_2d<double> fitter;
                    fitter.set_min_fit_length(min_fit_length);
                    fitter.set_rms_error_tol(0.05f);
                    fitter.add_curve(poly[0]);

                    vcl_vector<vgl_line_segment_2d<double> > segs;
                    fitter.fit();
                    segs= fitter.get_line_segs();

                    vgl_polygon<double> fitted_poly(1);
                    fitted_poly.push_back(segs[0].point1());
                    fitted_poly.push_back(segs[0].point2());

                    // See if segs intersects any bnd elements
                    for ( unsigned int i=1; i < segs.size() ; ++i)
                    {
                        fitted_poly.push_back(segs[i].point2());
                    }

                    closed_polys.push_back(fitted_poly);
                }
            }
        }


    }

    // vcl_map< double, vsol_spatial_object_2d_sptr >::iterator it;
    // for ( it = low_prob_curves.begin() ; it != low_prob_curves.end() ; ++it)
    // {
    //     if ( (*it).first < .13 )
    //     {
    //         input_vsol->remove_object((*it).second);
    //     }
    // }

    if ( remove_closed )
    {
        // Now erase closed contours removed
        for ( unsigned int c=0; c < objects_to_erase.size(); ++c)
        {
            input_vsol->remove_object(objects_to_erase[c]);

        }
    }


    /*********************** Shock Compute **********************************/
    // Grab output from shock computation
    vcl_vector<bpro1_storage_sptr> shock_results;
    {
        // 3) Create shock pro process and assign inputs 
        dbsk2d_compute_ishock_process shock_pro;

        shock_pro.clear_input();
        shock_pro.clear_output();

        // Default
        shock_pro.parameters()->set_value( "-rms" , (float)0.1 );

        shock_pro.add_input(frame_image);
        shock_pro.add_input(input_vsol);

        // Set params
        status = shock_pro.execute();
        shock_pro.finish();

        // If ishock status is bad we will keep iterating with noise 
        // till we get a valid shock computation otherwise call it quits
        if (!status)
        {
            // Add noise to parameter set
            shock_pro.parameters()->set_value("-b_noise",true);

            // Clean up before we start running
            shock_pro.clear_input();
            shock_pro.clear_output();

            unsigned int i(0);
            unsigned int num_iterations = 5;

            for ( ; i < num_iterations; ++i)
            {
                vcl_cout<<vcl_endl;
                vcl_cout<<"************ Retry Compute Shock,iter: "
                        <<i+1<<" *************"<<vcl_endl;

                // Add inputs
                shock_pro.add_input(frame_image);
                shock_pro.add_input(input_vsol);

                // Kick off process again
                status = shock_pro.execute();
                shock_pro.finish();

                if ( status )
                {
                    // We have produced valid shocks lets quit
                    break;

                }

                // Clean up after ourselves
                shock_pro.clear_input();
                shock_pro.clear_output();

            }
        }

        if ( status )
        {
            shock_results = shock_pro.get_output();

            // Clean up after ourselves
            shock_pro.clear_input();
            shock_pro.clear_output();

        }


    }

    vcl_string filename = output_folder+"/" + output_prefix +"_cgraph.dot";

    vcl_string stats_filename = output_folder+"/" + output_prefix +
                                "_stats_file.txt";
    vcl_ofstream stats_filestream(stats_filename.c_str());
    stats_filestream<<"Number of Contours: "<<vsol_list.size()<<vcl_endl;



    dbsk2d_shock_storage_sptr shock_storage;
    shock_storage.vertical_cast(shock_results[0]);


    // create_dependency_graph(shock_storage->get_ishock_graph(),
    //                         output_prefix);

    // Debug cost
    // debug_cost(shock_storage->get_ishock_graph(),
    //            output_prefix);

    // Debug frags
	vcl_string output_debug_frags = output_folder + "/" + output_prefix;
     debug_frags(shock_storage->get_ishock_graph(),
                 output_debug_frags);

    // pre_process_gap4(shock_storage);

    // Out of here

	dbsk2d_ishock_graph_sptr ishock_graph = shock_storage->get_ishock_graph();

	int count = 0;
	vcl_vector< dbsk2d_ishock_belm* > belms = ishock_graph->boundary()->belm_list();
	for (vcl_vector< dbsk2d_ishock_belm* >::iterator it = belms.begin(); it != belms.end(); it++)
	{
		if((*it)->type() == 0)
		{
			count++;
		}
	}

	stats_filestream << "Number of polyline nodes: " << count << vcl_endl;
	stats_filestream << "Number of ishock nodes: " << ishock_graph->all_nodes().size() << vcl_endl;
	stats_filestream<<vcl_endl;

	ishock_graph->ob_shocks();

    pre_process_contours(ishock_graph,
                         preprocess_threshold,
                         gap_distance,
                         remove_closed,
                         train);

    dbsk2d_containment_graph cgraph(ishock_graph,
                                    path_threshold,
                                    loop_cost,
                                    outside,
                                    train,
                                    debug,
                                    show_shock,
                                    quad);

    cgraph.construct_graph();
    cgraph.write_stats(stats_filestream);

    if ( debug )
    {
        cgraph.write_graph(filename);
    }


    if ( remove_closed )
    {
        for ( unsigned int c=0; c < closed_polys.size() ; ++c)
        {
            dbsk2d_transform_manager::Instance().write_output_region
                    (closed_polys[c]);
            dbsk2d_transform_manager::Instance().write_output_polygon
                    (closed_polys[c]);

            if ( train )
            {
                dbsk2d_transform_manager::Instance().write_stats_closed
                        (closed_polys[c]);
            }
        }
    }

    double vox_time = t.real()/1000.0;
    t.mark();
	timer = clock() - timer;
    vcl_cout<<"************ Time taken: "<<vox_time<<" sec"<<vcl_endl;
    stats_filestream<<"************ Time taken: "<<vox_time<<" sec"<<vcl_endl;
	stats_filestream << "*************New Time taken: " << (float)timer / CLOCKS_PER_SEC << "sec" << vcl_endl;
    stats_filestream.close();
    return status;
}

bool dbsk2d_compute_containment_graph_process::finish()
{
    return true;
}

void dbsk2d_compute_containment_graph_process::
create_dependency_graph(dbsk2d_ishock_graph_sptr ishock_graph,
                        vcl_string filename)
{

    // Detect transforms from new regions
    vcl_vector<vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bpoint*> >
            gap_pairs;
    vcl_vector<vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bline*> >
            gap4_pairs;

    vcl_map<vcl_pair<int,int> , vgl_point_2d<double> > gaps_train;
    {
        dbsk2d_ishock_gap_detector detector(ishock_graph);
        detector.detect_gap1(gap_pairs);

        for ( int i =0; i < gap_pairs.size() ; ++i)
        {
            vcl_pair<int,int>  id(gap_pairs[i].first->id(),
                                  gap_pairs[i].second->id());
            vgl_point_2d<double> midpoint =
                    _midPointPoint(gap_pairs[i].first->pt(),
                                   gap_pairs[i].second->pt());

            gaps_train[id]=midpoint;
        }


    }

    {
        dbsk2d_ishock_gap_detector detector(ishock_graph);
        detector.detect_gap4(gap4_pairs);


        for ( int i =0; i < gap4_pairs.size() ; ++i)
        {

            vgl_point_2d<double> f1=gap4_pairs[i].first->pt();
            vgl_point_2d<double> p1=gap4_pairs[i].second->s_pt()->pt();
            vgl_point_2d<double> p2=gap4_pairs[i].second->e_pt()->pt();

            vcl_pair<int,int> id;
            vgl_point_2d<double> midpoint;
            id.first = gap4_pairs[i].first->id();
            if ( vgl_distance(f1,p1) < vgl_distance(f1,p2) )
            {
                id.second = gap4_pairs[i].second->s_pt()->id();
                midpoint = _midPointPoint(gap4_pairs[i].first->pt(),
                                          gap4_pairs[i].second->s_pt()->pt());

            }
            else
            {
                id.second = gap4_pairs[i].second->e_pt()->id();
                midpoint = _midPointPoint(gap4_pairs[i].first->pt(),
                                          gap4_pairs[i].second->e_pt()->pt());

            }

            gaps_train[id]=midpoint;
        }
    }

    //Keep track of looops
    vcl_map<vcl_pair<int,int> , vgl_point_2d<double> > loops_train;

    vcl_set<vcl_pair<int,int> > loops_key;
    // Lets do loops first
    {

        vcl_vector<dbsk2d_ishock_belm*> belm_list = ishock_graph->
                boundary()->belm_list();

        vcl_vector<vsol_spatial_object_2d_sptr> vsol_list;

        for (unsigned int i=0;i < belm_list.size() ; ++i)
        {
            if ( belm_list[i]->is_a_point() )
            {
                dbsk2d_ishock_bpoint* bpoint =
                        dynamic_cast<dbsk2d_ishock_bpoint*>
                        (belm_list[i]);

                if ( bpoint->is_an_end_point() && bpoint->is_a_GUIelm())
                {

                    dbsk2d_ishock_loop_transform transformer(ishock_graph,
                                                             bpoint);

                    vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bpoint*>
                            pair=transformer.get_contour_pair();

                    vcl_pair<int,int> p1(pair.first->id(),pair.second->id());
                    vcl_pair<int,int> p2(pair.second->id(),pair.first->id());

                    vgl_point_2d<double> midpoint;
                    if ( !(loops_key.count(p1) || loops_key.count(p2)))
                    {

                        vcl_vector< vgl_point_2d<double> > ex_pts=
                                transformer.get_ordered_contour();

                        vcl_vector<vsol_point_2d_sptr> pts;

                        for ( unsigned int i=0; i < ex_pts.size() ; ++i)
                        {
                            vsol_point_2d_sptr pt=new vsol_point_2d(ex_pts[i]);

                            pts.push_back(pt);
                        }

                        dbsol_interp_curve_2d_sptr c =
                                new dbsol_interp_curve_2d();
                        dbsol_curve_algs::interpolate_linear(c.ptr(),
                                                             pts, false);

                        vsol_point_2d_sptr pt=c->point_at(c->length()/2.0);

                        midpoint.set(pt->x(),pt->y());

                        loops_train[p1]=midpoint;
                    }

                    loops_key.insert(p1);
                    loops_key.insert(p2);




                }
                else if ( bpoint->nLinkedElms()>= 6 )
                {
                    belm_list::iterator curB = bpoint->LinkedBElmList.begin();
                    for(; curB!=bpoint->LinkedBElmList.end(); ++curB)
                    {
                        if ( (*curB)->is_a_line())
                        {
                            dbsk2d_ishock_bline* bline=
                                    dynamic_cast<dbsk2d_ishock_bline*>
                                    (*curB);

                            dbsk2d_ishock_loop_transform transformer(
                                    ishock_graph,
                                    bpoint,
                                    bline);

                            vcl_pair<dbsk2d_ishock_bpoint*,
                                    dbsk2d_ishock_bpoint*>
                                    pair=transformer.get_contour_pair();

                            vcl_pair<int,int> p1(pair.first->id(),
                                                 pair.second->id());
                            vcl_pair<int,int> p2(pair.second->id(),
                                                 pair.first->id());

                            vgl_point_2d<double> midpoint;

                            if ( !(loops_key.count(p1) ||
                                   loops_key.count(p2)))
                            {

                                vcl_vector< vgl_point_2d<double> > ex_pts=
                                        transformer.get_ordered_contour();

                                vcl_vector<vsol_point_2d_sptr> pts;

                                for ( unsigned int i=0; i < ex_pts.size() ; ++i)
                                {
                                    vsol_point_2d_sptr pt=new vsol_point_2d(
                                            ex_pts[i]);

                                    pts.push_back(pt);
                                }

                                dbsol_interp_curve_2d_sptr c =
                                        new dbsol_interp_curve_2d();
                                dbsol_curve_algs::interpolate_linear(c.ptr(),
                                                                     pts,
                                                                     false);

                                vsol_point_2d_sptr pt=c->point_at(c->length()/
                                                                  2.0);

                                midpoint.set(pt->x(),pt->y());
                                loops_train[p1]=midpoint;
                            }

                            loops_key.insert(p1);
                            loops_key.insert(p2);
                        }
                    }
                }
            }

        }
    }

    vcl_string out=filename + "_depend_graph.txt";
    vcl_ofstream stream(out.c_str());

    vcl_map<vcl_pair<int,int> , vgl_point_2d<double> >::iterator it;
    for ( it = gaps_train.begin() ; it != gaps_train.end() ; ++it)
    {
        vcl_pair<int,int> orig=(*it).first;
        vgl_point_2d<double> p1=(*it).second;
        vcl_map<vcl_pair<int,int> , vgl_point_2d<double> >::iterator bit;
        bit=it;
        ++bit;
        stream<<p1.x()<<" "<<p1.y();
        for ( ; bit != gaps_train.end() ; ++bit)
        {
            vcl_pair<int,int> compare=(*bit).first;
            vgl_point_2d<double> p2=(*bit).second;
            if ( orig.first == compare.first ||
                 orig.first == compare.second ||
                 orig.second == compare.first ||
                 orig.second == compare.second )
            {
                stream<<" "<<p2.x()<<" "<<p2.y();
            }

        }

        vcl_map<vcl_pair<int,int> , vgl_point_2d<double> >::iterator lit;
        for ( lit=loops_train.begin() ; lit != loops_train.end() ; ++lit)
        {
            vcl_pair<int,int> compare=(*lit).first;
            vgl_point_2d<double> p2=(*lit).second;
            if ( orig.first == compare.first ||
                 orig.first == compare.second ||
                 orig.second == compare.first ||
                 orig.second == compare.second )
            {
                stream<<" "<<p2.x()<<" "<<p2.y();
            }

        }

        stream<<vcl_endl;

    }

}

void dbsk2d_compute_containment_graph_process::
debug_cost(dbsk2d_ishock_graph_sptr ishock_graph,
           vcl_string filename)
{

    // Write out boundary
    {
        vcl_string bnd_file=filename +"_pre_process.bnd";
        dbsk2d_ishock_transform temp_trans(ishock_graph,
                                           dbsk2d_ishock_transform::LOOP);
        temp_trans.write_boundary(bnd_file);
    }

    // Lets do loops first
    {
        vcl_string loop_filename = filename + "_loops.txt";
        vcl_ofstream loop_stream(loop_filename.c_str());

        vcl_vector<dbsk2d_ishock_belm*> belm_list = ishock_graph->
                boundary()->belm_list();

        //Keep track of looops
        vcl_set<vcl_pair<int,int> > loops_train;

        vcl_vector<vsol_spatial_object_2d_sptr> vsol_list;

        for (unsigned int i=0;i < belm_list.size() ; ++i)
        {
            if ( belm_list[i]->is_a_point() )
            {
                dbsk2d_ishock_bpoint* bpoint =
                        dynamic_cast<dbsk2d_ishock_bpoint*>
                        (belm_list[i]);

                if ( bpoint->is_an_end_point() && bpoint->is_a_GUIelm())
                {

                    dbsk2d_ishock_loop_transform transformer(ishock_graph,
                                                             bpoint);

                    vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bpoint*>
                            pair=transformer.get_contour_pair();

                    vcl_pair<int,int> p1(pair.first->id(),pair.second->id());
                    vcl_pair<int,int> p2(pair.second->id(),pair.first->id());

                    if ( !(loops_train.count(p1) || loops_train.count(p2)))
                    {
                        loop_stream<<bpoint->pt().x()<<" "
                                   <<bpoint->pt().y()<<" "
                                   <<transformer.likelihood()
                                   <<vcl_endl;

                        // Add in contours for front 
                        vsol_spatial_object_2d_sptr obj=
                                new vsol_polyline_2d;

                        vsol_polyline_2d* curve=obj->cast_to_curve()
                                ->cast_to_polyline();

                        vcl_vector< vgl_point_2d<double> > ex_pts=
                                transformer.get_ordered_contour();

                        for ( unsigned int i=0; i < ex_pts.size() ; ++i)
                        {
                            vsol_point_2d_sptr pt=new vsol_point_2d(ex_pts[i]);

                            curve->add_vertex(pt);
                        }

                        vsol_list.push_back(obj);

                    }


                    loops_train.insert(p1);
                    loops_train.insert(p2);


                }
                else if ( bpoint->nLinkedElms()>= 6 )
                {
                    belm_list::iterator curB = bpoint->LinkedBElmList.begin();
                    for(; curB!=bpoint->LinkedBElmList.end(); ++curB)
                    {

                        if ( (*curB)->is_a_line())
                        {
                            dbsk2d_ishock_bline* bline=
                                    dynamic_cast<dbsk2d_ishock_bline*>
                                    (*curB);

                            dbsk2d_ishock_loop_transform transformer(
                                    ishock_graph,
                                    bpoint,
                                    bline);

                            vcl_pair<dbsk2d_ishock_bpoint*,
                                    dbsk2d_ishock_bpoint*>
                                    pair=transformer.get_contour_pair();

                            vcl_pair<int,int> p1(pair.first->id(),
                                                 pair.second->id());
                            vcl_pair<int,int> p2(pair.second->id(),
                                                 pair.first->id());

                            if ( !(loops_train.count(p1) ||
                                   loops_train.count(p2)))
                            {
                                loop_stream<<bpoint->pt().x()<<" "
                                           <<bpoint->pt().y()<<" "
                                           <<transformer.likelihood()
                                           <<vcl_endl;

                                vsol_spatial_object_2d_sptr obj=
                                        new vsol_polyline_2d;

                                vsol_polyline_2d* curve=obj->cast_to_curve()
                                        ->cast_to_polyline();

                                vcl_vector<vgl_point_2d<double> > ex_pts=
                                        transformer.get_ordered_contour();

                                for ( unsigned int i=0; i < ex_pts.size() ; ++i)
                                {
                                    vsol_point_2d_sptr pt=new vsol_point_2d(
                                            ex_pts[i]);

                                    curve->add_vertex(pt);
                                }

                                vsol_list.push_back(obj);

                            }

                            loops_train.insert(p1);
                            loops_train.insert(p2);


                        }
                    }
                }
            }

        }

        loop_stream.close();

        vcl_string train_file=filename + "_train.cem";
        dbsol_save_cem(vsol_list, train_file);

        // write out train cons

    }



    // Lets do gaps first
    {

        vcl_string gap_filename = filename + "_gaps.txt";
        vcl_ofstream gap_stream(gap_filename.c_str());

        dbsk2d_ishock_gap_detector detector(ishock_graph);
        vcl_vector<
                vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bpoint*> > gap_pairs;

        detector.detect_gap1(gap_pairs);

        for ( unsigned int i=0; i < gap_pairs.size(); ++i)
        {
            int contour_id=vcl_min(gap_pairs[i].first->get_contour_id(),
                                   gap_pairs[i].second->get_contour_id());

            dbsk2d_ishock_gap_transform trans(ishock_graph,gap_pairs[i],
                                              contour_id);
            vcl_vector<vgl_point_2d<double> > samples=
                    trans.get_es_samples();
            gap_stream<<samples.size()<<" "<<trans.likelihood()<<vcl_endl;

            for ( int i=0; i < samples.size() ; ++i)
            {
                gap_stream<<samples[i].x()<<" "<<samples[i].y()<<vcl_endl;

            }

        }

        gap_stream.close();
    }

    {

        vcl_string gap_filename = filename + "_gap4s.txt";
        vcl_ofstream gap_stream(gap_filename.c_str());

        dbsk2d_ishock_gap_detector detector(ishock_graph);
        vcl_vector<
                vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bline*> > gap4_pairs;
        detector.detect_gap4(gap4_pairs);
        vcl_map<double,int> gap4_con_ids;

        for ( unsigned int i=0; i < gap4_pairs.size(); ++i)
        {
            int contour_id=this->get_contour(gap4_pairs[i].first)->get_id();
            dbsk2d_ishock_bpoint* anchor_pt =
                    dbsk2d_transform_manager::Instance()
                            .get_anchor_pt(gap4_pairs[i]);

            if ( anchor_pt->is_an_end_point())
            {
                continue;
            }

            dbsk2d_ishock_gap4_transform trans
                    (ishock_graph,gap4_pairs[i],anchor_pt,
                     contour_id);

            gap_stream<<gap4_pairs[i].first->pt().x()<<" "
                      <<gap4_pairs[i].first->pt().y()<<" "
                      <<anchor_pt->pt().x()<<" "
                      <<anchor_pt->pt().y()<<" "
                      <<trans.likelihood()<<vcl_endl;
        }

        gap_stream.close();
    }




}

void dbsk2d_compute_containment_graph_process::
debug_frags(dbsk2d_ishock_graph_sptr ishock_graph,
            vcl_string filename)
{

    // Write out boundary

    vcl_string shock_file  = filename +"_shock.cem";
    vcl_string cem_file    = filename + "_con.bnd";
    vcl_string atomic_frags= filename + "_atomic_frags";
    vcl_string coarse_poly = filename+"_coarse_fragments.txt";
    vcl_string coarse_shock= filename+"_csg.cem";

    dbsk2d_ishock_transform temp_trans(ishock_graph,
                                       dbsk2d_ishock_transform::LOOP);
    temp_trans.write_shock_boundary(shock_file);
    temp_trans.write_boundary(cem_file);

    //instantiate an empty coarse shock graph at this time
    //this will be properly defined once the shock is pruned
    dbsk2d_shock_graph_sptr coarse_shock_graph = new dbsk2d_shock_graph();

    // Write out coarse shock graph 
    dbsk2d_prune_ishock prune(ishock_graph, coarse_shock_graph);
    prune.prune(1.0);
    prune.compile_coarse_shock_graph();


    dbsk2d_ishock_grouping_transform grouper(ishock_graph);
    grouper.grow_regions();

    dbsk2d_ishock_grouping_transform grouper2(ishock_graph);
    grouper2.grow_coarse_regions();

    vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_edge*> >
            frag_edges = grouper.get_region_nodes();

    vcl_vector<vgl_polygon<double> > polygons;
    vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_edge*> >::iterator it;
    for ( it = frag_edges.begin() ; it != frag_edges.end() ; ++it)
    {
        int c=0;
        vcl_vector<dbsk2d_ishock_edge*> edges=(*it).second;
        for ( int b=0; b < edges.size() ; ++b)
        {
            if ( edges[b]->isHidden())
            {
                ++c;
            }
        }

        if ( c == edges.size() )
        {
            continue;
        }
        // Grab polygon fragment
        vgl_polygon<double> poly;
        grouper.polygon_fragment((*it).first,poly);

        polygons.push_back(poly);
    }

    temp_trans.write_fragments(filename,polygons);


    vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_edge*> >
            frag_edges2 = grouper2.get_region_nodes();

    vcl_vector<vgl_polygon<double> > type_1_polygons;
    vcl_vector<vgl_polygon<double> > type_2_polygons;
    vcl_vector<vgl_polygon<double> > atomic_polygons;

    vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_edge*> >::iterator it2;
    for ( it2 = frag_edges2.begin() ; it2 != frag_edges2.end() ; ++it2)
    {

        vcl_vector<dbsk2d_ishock_edge*> vec=(*it2).second;

        bool flag=false;
        for ( int i=0; i < vec.size() ; ++i)
        {
            flag=grouper2.shock_from_endpoint(vec[0]);
            if ( flag )
            {
                break;
            }
        }


        // Grab polygon fragment
        vgl_polygon<double> poly2;
        grouper2.polygon_fragment((*it2).first,poly2);
        atomic_polygons.push_back(poly2);

        if ( !flag )
        {
            type_1_polygons.push_back(poly2);
        }
        else
        {

            type_2_polygons.push_back(poly2);
        }
    }

    vcl_string type1=filename+"_type1_frags.txt";
    vcl_string type2=filename+"_type2_frags.txt";
    vcl_string type3=filename+"_type3_frags.txt";
    vcl_string type4=filename+"_csg_nodes.txt";
    vcl_string type5=filename+"_type_edges.txt";

    temp_trans.write_polygons(type1,type_1_polygons);
    temp_trans.write_polygons(type2,type_2_polygons);
    temp_trans.write_polygons(type3,atomic_polygons);
    grouper2.write_out_file(type4);
    grouper2.write_out_edges(type5);

    vcl_ofstream file_stream(coarse_poly.c_str());

    vcl_vector<vsol_spatial_object_2d_sptr> vsol_list;

    //draw the edges fragments first
    for ( dbsk2d_shock_graph::edge_iterator curE =
            coarse_shock_graph->edges_begin();
          curE != coarse_shock_graph->edges_end();
          curE++ )
    {
        dbsk2d_shock_edge_sptr selm = (*curE);

        //return if no fragment has been computed 
        if (!selm->shock_fragment())
        {
            continue;
        }

        file_stream<<selm->shock_fragment()->ex_pts().size()<<vcl_endl;

        for( unsigned int i = 0 ; i < selm->shock_fragment()->ex_pts().size()
                ; i++ ) {
            file_stream<<selm->shock_fragment()->ex_pts()[i].x()<<" "<<
                       selm->shock_fragment()->ex_pts()[i].y() <<" ";
        }
        file_stream<<vcl_endl;

        // Add in contours for front 
        vsol_spatial_object_2d_sptr obj=
                new vsol_polyline_2d;

        vsol_polyline_2d* curve=obj->cast_to_curve()->cast_to_polyline();

        vcl_vector<vgl_point_2d<double> > ex_pts= selm->ex_pts();


        for ( unsigned int i=0; i < ex_pts.size() ; ++i)
        {
            vsol_point_2d_sptr pt=new vsol_point_2d(ex_pts[i]);

            curve->add_vertex(pt);
        }

        vsol_list.push_back(obj);
    }

    dbsol_save_cem(vsol_list, coarse_shock);

    //then draw the node fragments

    for ( dbsk2d_shock_graph::vertex_iterator curN =
            coarse_shock_graph->vertices_begin();
          curN != coarse_shock_graph->vertices_end();
          curN++ )
    {
        dbsk2d_shock_node_sptr snode = (*curN);

        //traverse the descriptor list and draw the shock fragments for the 
        //degenerate descriptors
        vcl_list<dbsk2d_shock_node_descriptor>::iterator p_itr =
                snode->descriptor_list().begin();
        for (; p_itr != snode->descriptor_list().end(); ++ p_itr)
        {
            dbsk2d_shock_node_descriptor cur_descriptor = (*p_itr);

            if (!cur_descriptor.fragment.ptr())
            {
                continue;
            }

            file_stream<<cur_descriptor.fragment->ex_pts().size()<<vcl_endl;

            for( unsigned int i = 0 ; i <
                                      cur_descriptor.fragment->ex_pts().size() ; i++ ) {
                file_stream<<cur_descriptor.fragment->ex_pts()[i].x()<<" "<<
                           cur_descriptor.fragment->ex_pts()[i].y()<<" ";
            }

            file_stream<<vcl_endl;

        }
    }
    file_stream.close();

}

void dbsk2d_compute_containment_graph_process::
pre_process_contours(dbsk2d_ishock_graph_sptr ishock_graph,
                     double preprocess_threshold,
                     double gap_distance,
                     bool remove_closed,
                     bool train)
{

    vcl_vector<dbsk2d_ishock_belm*> belm_list = ishock_graph->boundary()->belm_list();
//    vcl_ofstream outfp("11.txt", vcl_ios::out);
//    for (vcl_vector< dbsk2d_ishock_belm* >::iterator it = belm_list.begin(); it != belm_list.end(); it++)
//    {
//        if ((*it)->type() == 0)
//        {
//            outfp << (*it)->id() << " ";
//            outfp << (*it)->type() << " ";
//            vgl_point_2d<double> pt = ((dbsk2d_ishock_bpoint*)(*it))->pt();
//            outfp << pt.x() << " " << pt.y() << " ";
//            outfp << vcl_endl;
//        }
//    }
//    outfp.close();
    for (unsigned int i=0;i < belm_list.size() ; ++i)
    {
        if ( belm_list[i]->is_a_point() )
        {
            dbsk2d_ishock_bpoint* bpoint =
                    dynamic_cast<dbsk2d_ishock_bpoint*>
                    (belm_list[i]);

            if ( bpoint->is_an_end_point() && bpoint->is_a_GUIelm())
            {
//                vcl_cout << "In pre_process_contours: In bpoint->is_an_end_point()" << vcl_endl;
//                vcl_cout << "bpoint id: " << bpoint->id() << vcl_endl;
                dbsk2d_ishock_loop_transform transformer(ishock_graph,
                                                         bpoint);

                vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bpoint*> pair =
                        transformer.get_contour_pair();
                double length_of_gap=vgl_distance(pair.first->pt(),
                                                  pair.second->pt());

                if ( length_of_gap < 1.0 )
                {
                    transformer.execute_transform();
                }

            }
        }

    }

    dbsk2d_ishock_gap_detector detector(ishock_graph);
    vcl_vector<
            vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bpoint*> > gap_pairs;

    detector.detect_gap1(gap_pairs);

    vcl_map<double,dbsk2d_ishock_transform_sptr> transforms;
    vcl_map<int,unsigned int> gaps_visited;
    for ( unsigned int i=0; i < gap_pairs.size(); ++i)
    {
        int contour_id=vcl_min(gap_pairs[i].first->get_contour_id(),
                               gap_pairs[i].second->get_contour_id());

        dbsk2d_ishock_transform_sptr trans = new
                dbsk2d_ishock_gap_transform(ishock_graph,gap_pairs[i],contour_id);
        transforms[trans->likelihood()]=trans;

        if ( gaps_visited.count(gap_pairs[i].first->id()) == 0 )
        {
            gaps_visited[gap_pairs[i].first->id()]=1;
        }
        else
        {
            gaps_visited[gap_pairs[i].first->id()]=
                    gaps_visited[gap_pairs[i].first->id()]+1;;
        }

        if ( gaps_visited.count(gap_pairs[i].second->id()) == 0 )
        {
            gaps_visited[gap_pairs[i].second->id()]=1;
        }
        else
        {
            gaps_visited[gap_pairs[i].second->id()]=
                    gaps_visited[gap_pairs[i].second->id()]+1;;
        }
    }

    vcl_map<int,vcl_vector<dbsk2d_bnd_contour_sptr> > gap_trans;
    vcl_map<double,dbsk2d_ishock_transform_sptr>::reverse_iterator it;
    for ( it = transforms.rbegin() ; it != transforms.rend() ; ++it)
    {
        vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bpoint*> pair =
                (*it).second->get_contour_pair();

        double length_of_gap=vgl_distance(pair.first->pt(),pair.second->pt());
        if ( (pair.first->is_an_end_point() && pair.second->is_an_end_point())
             &&
             (pair.first->is_a_GUIelm() && pair.second->is_a_GUIelm()))
        {
            if ( this->get_contour(pair.first)->get_id() ==
                 this->get_contour(pair.second)->get_id())
            {
                if ( gaps_visited[pair.first->id()] > 1 ||
                     gaps_visited[pair.second->id()] > 1 )
                {
                    continue;
                }

                dbsk2d_ishock_gap_transform* gap_transform
                        = (dbsk2d_ishock_gap_transform*)((*it).second.ptr());
                vcl_vector<vgl_point_2d<double> > gap_filler
                        = gap_transform->get_es_samples();

                dbsk2d_ishock_loop_transform loop_trans(
                        ishock_graph,
                        pair.first);
                vcl_vector<vgl_point_2d<double> > ordered_contour=
                        loop_trans.get_ordered_contour();

                vcl_vector<vgl_point_2d<double> > final_contour;
                if ( gap_filler[0] == ordered_contour[0])
                {
                    vcl_vector<vgl_point_2d<double> >::reverse_iterator rit;
                    for ( rit = gap_filler.rbegin(); rit != gap_filler.rend();
                          ++rit)
                    {
                        final_contour.push_back((*rit));
                    }
                }
                else
                {
                    vcl_vector<vgl_point_2d<double> >::iterator fit;
                    for ( fit = gap_filler.begin(); fit != gap_filler.end();
                          ++fit)
                    {
                        final_contour.push_back((*fit));
                    }
                }

                for ( unsigned int s=0; s < ordered_contour.size() ; ++s)
                {
                    final_contour.push_back(ordered_contour[s]);

                }

                vgl_polygon<double> poly(final_contour);

                dbsk2d_bnd_contour_sptr con1 = this->get_contour(pair.first);
                if ( gap_trans.count(con1->get_id()))
                {

                    if ( remove_closed )
                    {
                        dbsk2d_transform_manager::Instance().write_output_region
                                (gap_trans[con1->get_id()],gap_filler);
                        dbsk2d_transform_manager::Instance().
                                write_output_polygon
                                (poly);
                        if ( train )
                        {
                            dbsk2d_transform_manager::Instance()
                                    .write_stats_closed
                                            (poly);
                        }
                    }

                }
                else
                {

                    if ( remove_closed )
                    {
                        vcl_vector<dbsk2d_bnd_contour_sptr> cons;
                        cons.push_back(con1);
                        dbsk2d_transform_manager::Instance().write_output_region
                                (cons,gap_filler);
                        dbsk2d_transform_manager::Instance()
                                .write_output_polygon
                                        (poly);
                        if ( train )
                        {
                            dbsk2d_transform_manager::Instance()
                                    .write_stats_closed
                                            (poly);
                        }
                    }


                }

                if ( remove_closed )
                {
                    loop_trans.execute_transform();
                }
                else
                {
                    (*it).second->execute_transform();
                }
            }
            else if ( (*it).first >= preprocess_threshold &&
                      length_of_gap <= gap_distance)
            {
                dbsk2d_bnd_contour_sptr con1 = this->get_contour(pair.first);
                dbsk2d_bnd_contour_sptr con2 = this->get_contour(pair.second);


                int contour_id(0);
                if ( gap_trans.count(con1->get_id()))
                {
                    contour_id=con1->get_id();
                }
                else if ( gap_trans.count(con2->get_id()))
                {
                    contour_id=con2->get_id();
                }
                else
                {
                    contour_id=vcl_min(con1->get_id(),con2->get_id());
                }

                (*it).second->execute_transform();

                dbsk2d_bnd_contour_sptr gap_con =
                        ishock_graph->boundary()->preproc_contours().back();

                if ( !gap_trans.count(contour_id))
                {
                    gap_trans[contour_id].push_back(con1);
                    gap_trans[contour_id].push_back(con2);
                    gap_trans[contour_id].push_back(gap_con);
                    gap_con->set_id(contour_id);
                    con1->set_id(contour_id);
                    con2->set_id(contour_id);
                }
                else
                {
                    gap_trans[contour_id].push_back(gap_con);
                    gap_con->set_id(contour_id);
                    if ( con1->get_id() == contour_id)
                    {
                        if ( gap_trans.count(con2->get_id()))
                        {
                            vcl_vector<dbsk2d_bnd_contour_sptr> vec =
                                    gap_trans[con2->get_id()];
                            for ( unsigned int i=0; i < vec.size() ; ++i)
                            {
                                vec[i]->set_id(contour_id);
                                gap_trans[contour_id].push_back(vec[i]);

                            }
                        }
                        else
                        {
                            gap_trans[contour_id].push_back(con2);
                            con2->set_id(contour_id);
                        }
                    }
                    else
                    {
                        if ( gap_trans.count(con1->get_id()))
                        {
                            vcl_vector<dbsk2d_bnd_contour_sptr> vec =
                                    gap_trans[con1->get_id()];
                            for ( unsigned int i=0; i < vec.size() ; ++i)
                            {
                                vec[i]->set_id(contour_id);
                                gap_trans[contour_id].push_back(vec[i]);

                            }
                        }
                        else
                        {
                            con1->set_id(contour_id);
                            gap_trans[contour_id].push_back(con1);
                        }
                    }
                }


            }
        }

    }

    // {
    //     dbsk2d_ishock_transform temp_trans(ishock_graph,
    //                                        dbsk2d_ishock_transform::GAP);
    //     temp_trans.write_boundary("pre_processed_contours_gap1.bnd");
    // }

    transforms.clear();

    vcl_vector<
            vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bline*> > gap4_pairs;
    detector.detect_gap4(gap4_pairs);
    vcl_map<double,int> gap4_con_ids;

    for ( unsigned int i=0; i < gap4_pairs.size(); ++i)
    {
        int contour_id=this->get_contour(gap4_pairs[i].first)->get_id();
        dbsk2d_ishock_bpoint* anchor_pt = dbsk2d_transform_manager::Instance()
                .get_anchor_pt(gap4_pairs[i]);

        if ( anchor_pt->is_an_end_point())
        {
            continue;
        }

        dbsk2d_ishock_transform_sptr trans = new
                dbsk2d_ishock_gap4_transform(ishock_graph,gap4_pairs[i],anchor_pt,
                                             contour_id);
        transforms[trans->likelihood()]=trans;
        gap4_con_ids[trans->likelihood()]=
                gap4_pairs[i].second->get_contour_id();
    }

    for ( it = transforms.rbegin() ; it != transforms.rend() ; ++it)
    {
        vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bpoint*> pair =
                (*it).second->get_contour_pair();
        dbsk2d_ishock_gap4_transform* trans=(dbsk2d_ishock_gap4_transform*)
                ((*it).second.ptr());
        double length_of_gap=trans->length_of_gap();
        if ( pair.first->is_an_end_point() )
        {
            if ( (*it).second->valid_transform() )
            {
                if ( (*it).first >= preprocess_threshold &&
                     length_of_gap <= gap_distance )
                {
                    bool flag=(*it).second->execute_transform();

                    if ( flag )
                    {
                        // dbsk2d_ishock_bpoint* point=
                        //     (*it).second->endpoint_in_elms(
                        //         gap4_con_ids[(*it).first]);

                        // if ( point && (point->id() != pair.second->id()) )
                        // {
                        //     dbsk2d_ishock_loop_transform loop_trans(
                        //         ishock_graph,
                        //         point);
                        //     loop_trans.execute_transform();

                        // }
                    }
                    else
                    {

                        (*it).second->recompute_full_shock_graph();
                    }
                }
            }
        }

    }

    // {
    //     dbsk2d_ishock_transform temp_trans(ishock_graph,
    //                                        dbsk2d_ishock_transform::GAP);
    //     temp_trans.write_boundary("pre_processed_contours_gap4.bnd");
    // }

    transforms.clear();
    gap_pairs.clear();

    // dbsk2d_ishock_grouping_transform grouper(ishock_graph);
    // grouper.grow_regions();

    // vcl_map<unsigned int, vcl_vector<dbsk2d_ishock_belm*> > 
    //     region_belms = grouper.get_region_belms();
    // vcl_map<unsigned int,vcl_set<int> >
    //     region_belms_ids=grouper.get_region_belms_ids();
    // vcl_map<unsigned int, vcl_vector<dbsk2d_ishock_belm*> > 
    //     degree_three_nodes = grouper.get_degree_three_nodes();
    // vcl_map<unsigned int,vcl_set<int> >
    //     degree_three_node_ids=grouper.get_degree_three_node_ids();

    // vcl_map<vcl_pair<int,int> , vcl_pair<dbsk2d_ishock_bpoint*,
    //     dbsk2d_ishock_bline* > > degree_three_loops;

    // vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_belm*> >::iterator git;
    // for ( git = degree_three_nodes.begin() ; git != degree_three_nodes.end() ; 
    //       ++git)
    // {

    //     if ( !grouper.region_within_image((*git).first))
    //     {
    //         continue;
    //     }

    //     dbsk2d_ishock_bpoint* bpoint(0);
    //     dbsk2d_ishock_bline* bline(0);

    //     if ( (*git).second.size() >= 4 )
    //     {
    //         vcl_set<int> base_set=degree_three_node_ids[(*git).first];
    //         vcl_set<int> bline_base_set=region_belms_ids[(*git).first];

    //         vcl_map<unsigned int,vcl_vector<dbsk2d_ishock_belm*> >::iterator
    //             nit;
    //         nit=git;
    //         ++nit;
    //         for ( ; nit != degree_three_nodes.end() ; ++nit)
    //         {
    //             vcl_set<int> test_set=degree_three_node_ids[(*nit).first];
    //             vcl_set<int> bline_test_set=region_belms_ids[(*nit).first];

    //             vcl_set<int> intersection;
    //             vcl_insert_iterator<vcl_set<int> > 
    //                 inserter(intersection,intersection.begin());

    //             vcl_set_intersection(base_set.begin(),
    //                                  base_set.end(),
    //                                  test_set.begin(),
    //                                  test_set.end(),
    //                                  inserter);


    //             if ( intersection.size() >= 2 )
    //             {
    //                 vcl_set<int>::iterator sit;
    //                 for ( sit = intersection.begin() ; 
    //                       sit != intersection.end(); ++sit)
    //                 {
    //                     vcl_vector<dbsk2d_ishock_belm*> belms=(*git).second;
    //                     for ( unsigned int b=0; b < belms.size() ; ++b)
    //                     {
    //                         if ( belms[b]->id() == (*sit) )
    //                         {
    //                             bpoint = (dbsk2d_ishock_bpoint*)(belms[b]);
    //                         }

    //                     }
    //                 }


    //                 belm_list::iterator bit = bpoint->LinkedBElmList.begin();

    //                 for ( ; bit != bpoint->LinkedBElmList.end() ; ++bit)
    //                 {
    //                     bline=(dbsk2d_ishock_bline*)(*bit);

    //                     if ( bline_test_set.count(bline->twinLine()->id()) &&
    //                          bline_base_set.count(bline->id()))
    //                     {
    //                         break;
    //                     }
    //                 }


    //             }
    //         }
    //     }

    //     if ( bline && bpoint )
    //     {
    //         vcl_pair<int,int> key(bpoint->id(),bline->id());
    //         degree_three_loops[key]=vcl_make_pair(bpoint,bline);
    //     }
    // }

    // vcl_map<vcl_pair<int,int>,
    //     vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bline*> >::iterator lit;
    // for ( lit = degree_three_loops.begin() ; lit != degree_three_loops.end() ; 
    //       ++lit)
    // {

    //     dbsk2d_ishock_loop_transform transformer(ishock_graph,
    //                                              (*lit).second.first,
    //                                              (*lit).second.second);

    //     transformer.execute_transform();
    //     break;
    // }

    // {
    //     dbsk2d_ishock_transform temp_trans(ishock_graph,
    //                                        dbsk2d_ishock_transform::LOOP);
    //     temp_trans.write_boundary("pre_processed_contours_deg3_loop.bnd");
    // }


}

dbsk2d_bnd_contour_sptr dbsk2d_compute_containment_graph_process::
get_contour(dbsk2d_ishock_bpoint* bp)
{

    dbsk2d_bnd_vertex* vertex=bp->bnd_vertex();
    edge_list edges;
    vertex->edges(edges);

    edge_list::iterator it;
    if ( edges.size() == 2)
    {
        edges.erase(edges.begin());
    }
    it=edges.begin();

    const vcl_list< vtol_topology_object * > *
            superiors  = (*it)->superiors_list();
    vcl_list<vtol_topology_object*>::const_iterator tit;
    tit=(*superiors).begin();

    return (dbsk2d_bnd_contour*)(*tit);

}


void dbsk2d_compute_containment_graph_process::pre_process_gap4(
        dbsk2d_shock_storage_sptr output_shock)
{

    dbsk2d_ishock_graph_sptr ishock_graph=output_shock->get_ishock_graph();
    dbsk2d_boundary_sptr boundary=ishock_graph->boundary();

    vcl_vector<
            vcl_pair<dbsk2d_ishock_bpoint*,dbsk2d_ishock_bline*> > gap4_pairs;
    {
        dbsk2d_ishock_gap_detector detector(ishock_graph);
        detector.detect_gap4(gap4_pairs);
    }

    output_shock->set_ishock_graph(0);
    ishock_graph->remove_boundary();
    boundary->set_ishock_graph(0);
    ishock_graph->clear();
    ishock_graph=0;

    vil_image_resource_sptr img=dbsk2d_transform_manager::Instance().
            get_image();

    vcl_cout<<"Stats before splitting"<<vcl_endl;

    vcl_cout<<"Numb contours: "<<
            boundary->preproc_contours().size()+
            boundary->scratch_contours().size()
            <<vcl_endl;
    vcl_cout<<"Numb edges: "<<boundary->all_edges().size()
            <<vcl_endl;
    vcl_cout<<"Numb belms: "<<boundary->belm_list().size()
            <<vcl_endl;

    vcl_map<int, vcl_vector<dbsk2d_bnd_edge_sptr> > hash_map;
    for ( unsigned int i=0; i < gap4_pairs.size(); ++i)
    {

        dbsk2d_ishock_bpoint* bpt=gap4_pairs[i].first;
        dbsk2d_ishock_bline* bline=gap4_pairs[i].second;

        vgl_line_segment_2d<double> line(bline->s_pt()->pt(),
                                         bline->e_pt()->pt());

        vgl_point_2d<double> closest_pt = vgl_closest_point
                (line,bpt->pt());

        vgl_point_2d<double> pt1=bpt->pt();
        vgl_point_2d<double> pt2=closest_pt;


        if ( pt1.x() < 0 || pt1.y() < 0 ||
             pt1.x() >= img->ni() || pt1.y() >= img->nj() ||
             pt2.x() < 0 || pt2.y() < 0 ||
             pt2.x() >= img->ni() || pt2.y() >= img->nj())

        {
            continue;

        }

        dbsk2d_bnd_edge_sptr edge=0;

        if ( hash_map.count(bline->id()))
        {
            vcl_vector<dbsk2d_bnd_edge_sptr> edges=hash_map[bline->id()];
            for ( int i=0; i < edges.size() ; ++i)
            {

                if ( edges[i]->superiors_list()->size())
                {
                    vgl_line_segment_2d<double> test_line(
                            edges[i]->bnd_v1()->bpoint()->pt(),
                            edges[i]->bnd_v2()->bpoint()->pt());

                    if ( vgl_lineseg_test_point(closest_pt,
                                                test_line))
                    {
                        edge=edges[i];
                        break;
                    }
                }
            }

        }
        else if ( hash_map.count(bline->twinLine()->id()))
        {
            vcl_vector<dbsk2d_bnd_edge_sptr> edges=hash_map[
                    bline->twinLine()->id()];
            for ( int i=0; i < edges.size() ; ++i)
            {

                if ( edges[i]->superiors_list()->size())
                {
                    vgl_line_segment_2d<double> test_line(
                            edges[i]->bnd_v1()->bpoint()->pt(),
                            edges[i]->bnd_v2()->bpoint()->pt());

                    if ( vgl_lineseg_test_point(closest_pt,
                                                test_line))
                    {
                        edge=edges[i];
                        break;
                    }
                }
            }

        }
        else
        {
            edge=bline->bnd_edge();
        }


        if ( edge == 0 )
        {
            continue;
        }

        if ( vgl_distance(closest_pt,edge->bnd_v1()->bpoint()->pt()) < 0.1
             ||
             vgl_distance(closest_pt,edge->bnd_v2()->bpoint()->pt()) < 0.1
                )
        {
            continue;
        }


        dbsk2d_bnd_vertex_sptr anchor_pt=dbsk2d_bnd_utils::new_vertex(
                closest_pt,boundary);

        const vcl_list< vtol_topology_object * > *
                superiors  = edge->superiors_list();

        if ( superiors->size() == 0 )
        {
            break;
        }
        vcl_list<vtol_topology_object*>::const_iterator tit;
        tit=(*superiors).begin();

        dbsk2d_bnd_contour_sptr contour=(dbsk2d_bnd_contour*)(*tit);

        // now link all these vertices into a chain and save as a contour
        vcl_vector<dbsk2d_bnd_edge_sptr > bnd_edges;

        bnd_edges.push_back(dbsk2d_bnd_utils::new_line_between(
                edge->bnd_v1(),
                anchor_pt,
                boundary));

        bnd_edges.push_back(dbsk2d_bnd_utils::new_line_between(
                anchor_pt,
                edge->bnd_v2(),
                boundary));

        vcl_vector<signed char > directions(bnd_edges.size(), 1);

        bool output=contour->replace_edges(bnd_edges,
                                           directions,
                                           edge);

        bnd_edges[0]->left_bcurve()->set_GUIelm(true);
        bnd_edges[0]->right_bcurve()->set_GUIelm(true);

        bnd_edges[1]->left_bcurve()->set_GUIelm(true);
        bnd_edges[1]->right_bcurve()->set_GUIelm(true);

        hash_map[bline->id()].push_back(bnd_edges[0]);
        hash_map[bline->id()].push_back(bnd_edges[1]);

        if ( output == false)
        {
            break;
        }
    }

    boundary->update_belm_list();

    vcl_cout<<"Stats after splitting"<<vcl_endl;

    vcl_cout<<"Numb contours: "<<
            boundary->preproc_contours().size()+
            boundary->scratch_contours().size()
            <<vcl_endl;
    vcl_cout<<"Numb edges: "<<boundary->all_edges().size()
            <<vcl_endl;
    vcl_cout<<"Numb belms: "<<boundary->belm_list().size()
            <<vcl_endl;

    vcl_vector<vsol_spatial_object_2d_sptr> vsol_list;

    bnd_contour_list all_contours;
    boundary->all_contours(all_contours);
    unsigned int id=1;
    for (bnd_contour_list::iterator cit =
            all_contours.begin(); cit != all_contours.end(); ++cit)
    {

        // Add in contours for front 
        vsol_spatial_object_2d_sptr obj=
                new vsol_polyline_2d;

        vsol_polyline_2d* curve=obj->cast_to_curve()->cast_to_polyline();

        dbsk2d_bnd_contour_sptr con=(*cit);
        for ( int e=0; e <= con->num_edges() ; ++e)
        {
            vsol_point_2d_sptr pt=new vsol_point_2d(
                    con->bnd_vertex(e)->bpoint()->pt());
            curve->add_vertex(pt);
        }

        curve->set_id(id);
        vsol_list.push_back(curve);
        ++id;
    }

    // Now destroy boundary
    output_shock->set_boundary(0);
    output_shock->set_image(0);
    boundary=0;

    /*********************** Shock Compute **********************************/
    // Grab output from shock computation

    while ( true )
    {
        dbsk2d_boundary_sptr new_boundary=
                dbsk2d_create_boundary (
                        vsol_list,      // vsol objects
                        false,          // bool override_default_partitioning,
                        0,0,            // xmin,ymin
                        1,1,            // int num_rows, int num_cols,
                        1000.0f,1000.0f,// float cell_width, float cell_height,
                        true,    //bool preprocess_boundary,
                        true);    //bool break_long_lines,

        dbsk2d_ishock_graph_sptr ishock_new_graph
                =dbsk2d_compute_ishocks(new_boundary);

        if ( ishock_new_graph )
        {
            output_shock->set_ishock_graph(ishock_new_graph);
            output_shock->set_boundary(new_boundary);

            break;
        }
        else
        {
            // add noise to vsol list

            output_shock->set_boundary(0);
            output_shock->set_ishock_graph(0);

            vcl_cout<<"After preprocessing adding noise"<<vcl_endl;

            for (unsigned int v=0; v < vsol_list.size() ; ++v)
            {
                vsol_polyline_2d_sptr curve=vsol_list[v]->cast_to_curve()->
                        cast_to_polyline();


                vnl_random mz_random;
                mz_random.reseed((unsigned long)time(NULL));
                float noise_radius=0.002f;

                if ( curve->vertex(0)->get_p() == curve->vertex(
                        curve->size()-1)->get_p())
                {
                    // closed curve ignore
                    continue;
                }

                for ( int c=0; c < curve->size() ; ++c)
                {
                    vsol_point_2d_sptr vertex=curve->vertex(c);

                    vgl_point_2d<double> point=vertex->get_p();

                    double x=point.x();
                    double y=point.y();
                    double rand_x = mz_random.drand32(1.0);
                    x += 2.0*noise_radius*(rand_x-0.5);
                    double rand_y = mz_random.drand32(1.0);
                    y += 2.0*noise_radius*(rand_y-0.5);
                    vertex->set_x(x);
                    vertex->set_y(y);
                }
            }

        }
    }
    // dbsk2d_ishock_transform transform(output_shock->get_ishock_graph(),
    //                                   dbsk2d_ishock_transform::LOOP);
    // {
    //     transform.write_shock_boundary("inserted_gap4s_shocks.cem");
    //     transform.write_boundary("inserted_gap4s_contours.bnd");
    // }

}
