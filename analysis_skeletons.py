
'''
General skeleton analysis
'''

import numpy as np

from fil_finder.pixel_ident import pix_identify, extremum_pts, \
    make_final_skeletons, recombine_skeletons, isolateregions
from fil_finder.length import init_lengths, pre_graph, \
    longest_path, prune_graph, main_length


def analyze_skeletons(skeleton, relintens_thresh=0.2, imgscale=1.,
                      skel_thresh=None, branch_thresh=None,
                      verbose=False, save_png=False, save_name=None):
    '''

    This function wraps most of the skeleton analysis. Several steps are
    completed here:
    *   isolatefilaments is run to separate each skeleton into its own
        array. If the skeletons are under the threshold set by
        size_thresh, the region is removed. An updated mask is
        also returned.
    *   pix_identify classifies each of the pixels in a skeleton as a
        body, end, or intersection point. See the documentation on find_filpix
        for a complete explanation. The function labels the branches and
        intersections of each skeletons.
    *   init_lengths finds the length of each branch in each skeleton and
        also returns the coordinates of each of these branches for use in
        the graph representation.
    *   pre_graph turns the skeleton structures into a graphing format
        compatible with networkx. Hubs in the graph are the intersections
        and end points, labeled as letters and numbers respectively. Edges
        define the connectivity of the hubs and they are weighted by their
        length.
    *   longest_path utilizes networkx.shortest_path_length to find the
        overall length of each of the filaments. The returned path is the
        longest path through the skeleton. If loops exist in the skeleton,
        the longest path is chosen (this shortest path algorithm fails when
        used on loops).
    *   extremum_pts returns the locations of the longest path's extent
        so its performance can be evaluated.
    *   final_lengths takes the path returned from longest_path and
        calculates the overall length of the filament. This step also acts
        as to prune the skeletons.
    *   final_analysis combines the outputs and returns the results for
        further analysis.

    Parameters
    ----------
    verbose : bool, optional
        Enables plotting.
    relintens_thresh : float, optional
        Relative intensity threshold for pruning. Sets the importance
        a branch must have in intensity relative to all other branches
        in the skeleton. Must be between (0.0, 1.0].
    skel_thresh : float, optional
        Manually set the minimum skeleton threshold. Overrides all
        previous settings.
    branch_thresh : float, optional
        Manually set the minimum branch length threshold. Overrides all
        previous settings.
    save_png : bool, optional
        Saves the plot made in verbose mode. Disabled by default.

    Attributes
    ----------
    filament_arrays : list of numpy.ndarray
        Contains individual arrays of each skeleton
    number_of_filaments : int
        The number of individual filaments.
    array_offsets : list
        A list of coordinates for each filament array.This will
        be used to recombine the final skeletons into one array.
    filament_extents : list
        This contains the coordinates of the initial and final
        position of the skeleton's extent. It may be used to
        test the performance of the shortest path algorithm.
    lengths : list
        Contains the overall lengths of the skeletons
    labeled_fil_arrays : list of numpy.ndarray
        Contains the final labeled versions of the skeletons.
    branch_properties : dict
        The significant branches of the skeletons have their length
        and number of branches in each skeleton stored here.
        The keys are: *filament_branches*, *branch_lengths*

    '''

    if relintens_thresh > 1.0 or relintens_thresh <= 0.0:
        raise ValueError(
            "relintens_thresh must be set between (0.0, 1.0].")

    # Minimum skeleton length
    if skel_thresh is None:
        skel_thresh = 0

    # Minimum branch length
    if branch_thresh is None:
        branch_thresh = 0

    isolated_filaments, num, offsets = \
        isolateregions(skeleton, size_threshold=skel_thresh,
                       pad_size=2, fill_hole=False)

    interpts, hubs, ends, filbranches, labeled_fil_arrays =  \
        pix_identify(isolated_filaments, num)

    branch_properties = \
        init_lengths(labeled_fil_arrays, filbranches,
                     offsets, np.ones_like(skeleton))

    # Add the number of branches onto the dictionary
    branch_properties["number"] = filbranches

    edge_list, nodes = \
        pre_graph(labeled_fil_arrays, branch_properties, interpts, ends)

    max_path, extremum, G = \
        longest_path(edge_list, nodes,
                     verbose=verbose,
                     save_png=save_png,
                     save_name=save_name,
                     skeleton_arrays=labeled_fil_arrays,
                     lengths=branch_properties["length"])

    updated_lists = \
        prune_graph(G, nodes, edge_list, max_path, labeled_fil_arrays,
                    branch_properties, branch_thresh,
                    relintens_thresh=relintens_thresh)

    labeled_fil_arrays, edge_list, nodes, branch_properties = \
        updated_lists

    filament_extents = \
        extremum_pts(labeled_fil_arrays, extremum, ends)

    length_output = main_length(max_path, edge_list, labeled_fil_arrays,
                                interpts,
                                branch_properties["length"], imgscale,
                                verbose=False, save_png=save_png,
                                save_name=save_name)

    lengths, filament_long_paths = length_output
    # Convert lengths to numpy array
    lengths = np.asarray(lengths)

    filament_final_skeletons =\
        make_final_skeletons(labeled_fil_arrays, interpts,
                             verbose=verbose, save_png=save_png,
                             save_name=save_name)

    # Convert branch lengths physical units
    for n in range(num):
        b_lengths = branch_properties["length"][n]
        branch_properties["length"][n] = \
            [imgscale * length for length in b_lengths]

    skeleton = \
        recombine_skeletons(filament_final_skeletons,
                            offsets, skeleton.shape,
                            1, verbose=True)

    skeleton_longpath = \
        recombine_skeletons(filament_long_paths,
                            offsets, skeleton.shape,
                            1, verbose=True)

    return lengths, skeleton, skeleton_longpath, filament_extents