#include <gurobi_c++.h>
#include<cmath>
#include<map>
#include <iterator>
#include <utility>
#include<set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <iostream>  //std::right, std::out, std::endl
#include <iomanip> //std::setw(int w)
#include <fstream>   // ifstream, ofstream
#include <sstream>   // istringstream, ostringstream
#include <vector>
//#include "get_input.h"


using namespace std;

struct my_node_property {
//    double un-normalized_weight;
    //std::vector<double> contrib_h;
    //int id;
//    double flow;

    // This label represents node index in kf format (first argument in the (*?*:*) format).
    int label;
    int path_count;
    int parent_branch_instruction;
};

struct my_edge_property {
//    double weight;
    int zero_one;

    // The following set for the current case of NZDDs having one edge on each edge will reduce to a single member set.
    std::set<int> labels; // labels of edge
};

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
my_node_property, //properties of nodes
my_edge_property, //properties of edges
boost::no_property //properties of graph
> Graph;
typedef std::pair<int, int> Edge;
typedef boost::graph_traits<Graph>::vertex_descriptor Node;

std::map<int, Graph::vertex_descriptor> V_global;
//std::map<int, Graph::vertex_descriptor> V;


map<int, vector<int>> get_input(std::string filename, std::vector<std::vector<int>>& instances, std::vector<int> & labels,  int* number_of_samples, int* number_of_base_hypotheses);
void graph_recursive_construct(Graph& g, vector<Graph::vertex_descriptor>& vertices,
        vector<Graph::vertex_descriptor>& v_sink_one, vector<Graph::vertex_descriptor>& v_sink_zero,
        map<int, vector<int>> branch_instructions, int branch_instruction_number, int parent_branch_instruction_number,
        int zero_one);

int main(int argc, char *argv[]){


    int number_of_samples, number_of_base_hypotheses;
    int i, j, k;
    double eps =0.0001;
    std::vector<std::vector<int>> instances;
    std::vector<int> labels;
    std::vector<int>::iterator itr;
    //double xvals[];

//    std::string filename = argv[1];
//    std::string filename = "rofk-100-p.kf";
    std::string filename = "small.kf";



    //get input and store constraints in C_pos and C_neg
    map<int, vector<int>> branch_instructions;
    branch_instructions = get_input(filename, instances, labels, &number_of_samples, &number_of_base_hypotheses);
    cout << "Size of branch_instructions:" << branch_instructions.size() << endl;
    cout <<"number_of_samples==" << number_of_samples << endl;
    cout <<"number_of_base_hypotheses=" << number_of_base_hypotheses << endl;
    //cout <<"|X_pos|=" << X_pos.size() << endl;
    //cout <<"|X_neg|=" << X_neg.size() << endl;

    //SECTION: constructing the graph

    Graph g;
    vector<Graph::vertex_descriptor> V, V_sink_one_edge, V_sink_zero_edge; // those that are to be connected to sink via
    // one edge and zero edge should be kept separate.
    auto iter = branch_instructions.end();
    iter--;
    int last_instruction_number = (iter->first);
    cout << "###### Last Instruction number:\t" << (iter->first) << endl;

    // adding sink node
    Graph::vertex_descriptor one_sink_node = add_vertex(g);
    g[one_sink_node].label = 0;
    g[one_sink_node].parent_branch_instruction = 1; // i.e. it's the sink node and branch instruction for it is certainly 1!

//    std::pair<std::map<std::string, int>::iterator, bool > result;
    auto result = V_global.insert(std::pair<int, Graph::vertex_descriptor>(1, one_sink_node));


    // recursively constructing the graph
    graph_recursive_construct(g, V,V_sink_one_edge, V_sink_zero_edge, branch_instructions, last_instruction_number, 0, 0); // Fourth element is the pointer to the desired

    // when the above command is executed, graph is constructed along with its one and zero edges. Just edges to the
    // sink node remain to be added. next step is to do the topological sorting. BUT I don't think that the sorting is
    // required, since when kf format file is read from the end, the resulting recursive graph constructed, is itself
    // sorted.
    //NOTE: So we only need to add one-edge and zero-edge of the nodes in the V_sink_one_edge and V_sink_zero_edge
    // vectors. Then transform the graph to the NZDD form with a single label on each edge.

    // adding the sink node
//    auto sink_node = add_vertex(g);
    // Adding the one-edges
//    for (auto first = V_sink_one_edge.begin(), last = V_sink_one_edge.end(); first != last; ++first) {
//        cout << "ONE-EDGE to the sink node was added. From " << (g[*first]).label << endl;
//
//        bool inserted = false;
//        Graph::edge_descriptor e;
//
//        boost::tie(e, inserted) = add_edge(*first, sink_node, g);
//        g[e].zero_one = 1;
//    }
//    for (auto first = V_sink_zero_edge.begin(), last = V_sink_zero_edge.end(); first != last; ++first) {
//        cout << "ZERO-EDGE to the sink node was added. From " << (g[*first]).label << endl;
//
//        bool inserted = false;
//        Graph::edge_descriptor e;
//
//        boost::tie(e, inserted) = add_edge(*first, sink_node, g);
//        g[e].zero_one = 0;
//    }

    //SECTION: transformation of the graph to NZDD

    // In this part we have to iterate over all edges and assign their parent's label to them.
    //NOTE: Only one-edges, get their parent's label and 0 (or empty set) will be assigned to zero-edges.

    //SECTION: calculating the m_e for each edge using dynamic programming
    std::pair<boost::adjacency_list<>::vertex_iterator,
    boost::adjacency_list<>::vertex_iterator> vs = boost::vertices(g);

    for (boost::adjacency_list<>::vertex_iterator v_itr = vs.first; v_itr != vs.second; ++v_itr){
        cout << "m_e for node #" << g[(*v_itr)].label << " from parent branch instruction #" <<
        g[(*v_itr)].parent_branch_instruction << " is: " << g[(*v_itr)].path_count << endl;
    }

    cout << "Program Finished!" << endl;

    return 0;

}

map<int, vector<int>> get_input(std::string filename, std::vector<vector<int>>& instances, std::vector<int> & labels,
        int* number_of_samples, int* number_of_base_hypotheses){


    //input stream
    ifstream reading_file;
    reading_file.open(filename, std::ios::in);
    cout <<"reading done" << std::endl;
    string line;

    vector<string> splitted;
    vector<string> term;
    vector<string>::iterator itr;
    map<int, vector<int>> branch_instructions;

    // should be cleared after each loop.
    vector<int> tmp;

    string delim("(~,:,?,)");
    int index_of_sample;
    int index_of_hypothesis;
    int num_line =0;
    int dimensionality=0;
//    double val;
    int val;

    // Note: Each entry in the kf format, represents an internal node of the ZDD.
    // Suppose we have 4: (~8?0:1), this means that branch Instruction 4 corresponds to node 8 and if in a sample
    // or equivalently path, we have this node (or base hypothesis, since nodes represent the elements of our family set)
    // we will move to branch instruction 0 and otherwise to branch instruction 1.

    // Very Important! In the above example, 4,0 and 1 are all branch instructions and can be interpreted as labels of
    // edges between nodes. These nodes' values simply represent the index of hypothesis!

    // If it's not clear to you, first read the README of zcomp + p. 205-206 from the Art of Programming 4a. Then read again
    // and you'll understand what exactly is the notion of label and node value.

    // Representation acquired, can be easily converted to an NZDD. Transferring the node label to its 1 edge and labeling
    // zero edges by empty set. Then we're ready to calculate the m_e s.
    while (std::getline(reading_file, line)) {
        boost::trim(line); // remove spaces of both edges
        boost::algorithm::split(splitted, line,  boost::is_space()); //split the line with a space
        labels.push_back(stoi(splitted[0])); // the first part corresponds to label

        int instruction_index = stoi(splitted[0]);

        cout << "Label= " << stoi(splitted[0]) << endl;
        for(itr = splitted.begin(), ++itr; itr != splitted.end(); itr++) {
            cout << "Next Line!" << endl;
            std::string st(*itr);
            boost::trim (st);
            cout << "st= " << st <<endl;
            boost::algorithm::split(term, st, boost::is_any_of(delim));
//            cout << "term= " << term[0] << endl;
//            cout << "term= " << term[1] << endl;
//            cout << "term= " << term[2] << endl;
//            cout << "term= " << term[3] << endl;
//            cout << "term= " << term[4] << endl;
            ///////////////////////////////////////////// Amin
            tmp.push_back(stoi(term[2]));
//            cout << "tmp:" << tmp.at(0) << endl;
            tmp.push_back(stoi(term[3]));
//            cout << "tmp:" << tmp.at(1) << endl;
            tmp.push_back(stoi(term[4]));
//            cout << "tmp:" << tmp.at(2) << endl;
            cout << "Hypothesis:" << term[2] << "\t one-edge:" << term[3] << "\t zero-edge:" << term[4] << endl;
            branch_instructions[instruction_index] = tmp;
            tmp.clear();
            ///////////////////////////////////////////// Amin
            index_of_hypothesis = stoi(term[2]);
            val = stoi(term[3]);
            //cout << "index_h="<<index_h << "val=" << val <<endl;
            //std::istringstream index_val_pair(term);
            //index_val_pair >> index_h >> ":" >> val ;
            dimensionality = max(dimensionality, index_of_hypothesis);
        }
        num_line++;
    }

    //initialization (rewinding)
    reading_file.clear();
    reading_file.seekg(0, std::ios::beg);

    instances.resize(num_line);

    index_of_sample = 0;
    while (std::getline(reading_file, line)) {
        boost::trim(line);
        //cout << line << endl;
        boost::split(splitted, line, boost::is_space()); //split the line with a space
        for(itr = splitted.begin(), ++itr; itr != splitted.end(); ++itr) {
            cout << "Second for loop!" << endl;
            std::string st(*itr);
            boost::trim (st);
            cout << "st= " << st <<endl;
            boost::split(term, st, boost::is_any_of(delim));
            index_of_hypothesis =stoi(term[2])-1;
            val = stoi(term[3]);
            //cout << *itr <<endl;
            //std::istringstream index_val_pair(*itr);
            //index_val_pair >> index_h >> ":" >> val ;
            instances[index_of_sample].push_back(index_of_hypothesis);
        }
        index_of_sample++;
    }


    *number_of_samples = index_of_sample;
    *number_of_base_hypotheses = dimensionality;

    return branch_instructions;
}

void graph_recursive_construct(Graph& g, vector<Graph::vertex_descriptor>& vertices_descriptors,
        vector<Graph::vertex_descriptor>& v_sink_one, vector<Graph::vertex_descriptor>& v_sink_zero,
        map<int, vector<int>> branch_instructions, int branch_instruction_number, int parent_branch_instruction_number,
        int zero_one){

    cout << "Node index: " << (branch_instructions[branch_instruction_number]).at(0) << "\tZero-edge: " <<
    (branch_instructions[branch_instruction_number]).at(1) << "\tOne-edge: " <<
    (branch_instructions[branch_instruction_number]).at(2) << endl;

    ///////////////////////////////// Recursive

    //NOTE: if some node is observed, then we should NOT add a vertex for it. We should just connect it to the previously
    // observed node and terminate recursion through this branch since if that node is observed previously, certainly
    // recursion has been called on its children in another recursion path. For this purpose we should keep a GLOBAL map
    // of observed nodes available to all instances of recursion calls being undertaken at each moment.
    // map<branch_instruction_number, vertex_descriptor> or
    // map<is_observed, branch_instruction_number, vertex_descriptor> could be used.

    //NOTE: To calculate the m_e, at the instant of creation of a vertex, a number reflecting the number of paths from
    // the root to this point is inherited from the parent. If at any future instance of time, an observed node is the
    // child of a another node, number of that node is added to that of this node's hence it gets updated. SO WE NEED TO
    // PASS THE PARENT'S BRANCH INSTRUCTION AS AN ARGUMENT TO THE RECURSIVE FUNCTION CALL.

    // if the node has not been added yet
    if(V_global.find(branch_instruction_number) == V_global.end()){ // not observed before

        Graph ::vertex_descriptor temporary_descriptor;
        temporary_descriptor = add_vertex(g);
        vertices_descriptors.push_back(temporary_descriptor);
        g[temporary_descriptor].label = (branch_instructions[branch_instruction_number]).at(0);
        g[temporary_descriptor].parent_branch_instruction = branch_instruction_number; // if it's a root it will be zero,
        // otherwise is also correct

        V_global.insert(std::pair<int, Graph ::vertex_descriptor>(branch_instruction_number, temporary_descriptor));
        // Nodes label is set to hypothesis node and later will be given to its one-edge.
        cout << "Node " << g[temporary_descriptor].label << " for branch instruction #" << branch_instruction_number <<
        " was added." << endl;

        // connecting this node to its parent and inheriting the path count from the parent
        if (parent_branch_instruction_number == 0){// i.e. the node is a root itself and no edge to parent is required.
            // path count is initialized to 1.
            g[V_global[branch_instruction_number]].path_count = 1;

        } else{ // i.e. the node is not a root and needs to be connected to its parent via a one or zero edge and
            // inherit the path count of its parent.
            //NOTE: Only add and inherit the path_count from the parent if they are connected by a one-edge.
            //TODO: What to do for zero-edges? For now I'll add for both zero and one edges since I don't see any
            // difference

            bool inserted = false;
            Graph::edge_descriptor e;

            auto parent_descriptor = V_global[parent_branch_instruction_number];
            boost::tie(e, inserted) = add_edge(parent_descriptor, temporary_descriptor, g);
            g[e].zero_one = zero_one;

            g[V_global[branch_instruction_number]].path_count = g[V_global[parent_branch_instruction_number]].path_count;


        }

        // Defining children of this node. Since it had not been observed, its children have not been observed, too!
//        Graph ::vertex_descriptor temporary_descriptor_one_child, temporary_descriptor_zero_child;

        // 1. if the zero-child is zero, do nothing and leave recursion and the descendant.
        // 2. if the zero-child is one, connect the ZERO-edge to SINK node and finish recursion. For this reason we need a vector of descriptors that
        //      need to be connected to the sink node at the end of recursion.
        // 3. otherwise continue to its descendants and continue the recursion.
        if((branch_instructions[branch_instruction_number]).at(1) == 0 || (branch_instructions[branch_instruction_number]).at(1) == 1){
            // case 1
            if ((branch_instructions[branch_instruction_number]).at(1) == 0){
                cout << "zero-child is zero and recursion stops." << endl;

                //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
                // nothing should be happened according to these nodes.
            }
                // case 2
            else {
                //TODO::DONE Add to sink_descriptors
                v_sink_zero.push_back(temporary_descriptor);

                // Adding the zero-edge to the sink node
                cout << "ZERO-EDGE to the sink node was added. From branch instruction #" <<
                branch_instruction_number << " , node #" << g[V_global[branch_instruction_number]].label << endl;

                bool inserted = false;
                Graph::edge_descriptor e;

                boost::tie(e, inserted) = add_edge(temporary_descriptor, V_global[1], g);
                g[e].zero_one = 0;

                //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
                // nothing should be happened according to these nodes.
            }
            // if goes to the following else, it means that recursion must be continued for zero-edge.
        }
            // case 3
        else {

            // there are non-zero and non-one children, hence we go through the recursion to consider their branch instructions
            auto zero_child_instruction_number = branch_instructions[branch_instruction_number].at(1);

//            bool inserted = false;
//            Graph::edge_descriptor e;
//
//            boost::tie(e, inserted) = add_edge(temporary_descriptor, temporary_descriptor_zero_child, g);
//            std::cout<< "connecting node #" << g[temporary_descriptor].label << " to the node (ZERO-EDGE) #" << g[temporary_descriptor_zero_child].label << std::endl;
//            g[e].zero_one = 0;

            graph_recursive_construct(g, vertices_descriptors, v_sink_one, v_sink_zero, branch_instructions,
                    zero_child_instruction_number, branch_instruction_number, 0);
        }

        // 1. if the one-child is zero, do nothing and leave recursion and the descendant.
        // 2. if the one-child is one, connect the ONE-edge to SINK node and finish recursion. For this reason we need a vector of descriptors that
        //      need to be connected to the sink node at the end of recursion.
        // 3. otherwise continue to its descendants and continue the recursion.
        if((branch_instructions[branch_instruction_number]).at(2) == 0 || (branch_instructions[branch_instruction_number]).at(2) == 1){
            // case 1
            if ((branch_instructions[branch_instruction_number]).at(2) == 0){
                cout << "one-child is zero and recursion stops." << endl;
                //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
                // nothing should be happened according to these nodes.
            }
                // case 2
            else {
                //TODO::DONE Add to sink_descriptors
                v_sink_one.push_back(temporary_descriptor);

                // Adding the zero-edge to the sink node
                cout << "ONE-EDGE to the sink node was added. From branch instruction #" <<
                     branch_instruction_number << " , node #" << g[V_global[branch_instruction_number]].label << endl;

                bool inserted = false;
                Graph::edge_descriptor e;

                boost::tie(e, inserted) = add_edge(temporary_descriptor, V_global[1], g);
                g[e].zero_one = 1;

                //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
                // nothing should be happened according to these nodes.


            }
            // if goes to the following else, it means that recursion must be continued for one-edge.
        }
            // case 3
        else {
            /*TODO::DONE Pay attention to the shifts in kf data. some start from 2!
         * because of this, we use map instead of a pure vector.
         * */

//            temporary_descriptor_one_child = add_vertex(g);
//            vertices_descriptors.push_back(temporary_descriptor_one_child);

            // there are non-zero and non-one children, hence we go through the recursion to consider their branch instructions
            auto one_child_instruction_number = branch_instructions[branch_instruction_number].at(2);
            //            auto one_child_instruction_number = branch_instructions[branch_instruction_number].at(2);

            // This is the one-edge.
//            g[temporary_descriptor_one_child].label = branch_instructions[one_child_instruction_number].at(0);

//            bool inserted = false;
//            Graph::edge_descriptor e;

//            boost::tie(e, inserted) = add_edge(temporary_descriptor, temporary_descriptor_one_child, g);
//            std::cout<< "connecting node #" << g[temporary_descriptor].label << " to the node (ONE-EDGE) #" << g[temporary_descriptor_one_child].label << std::endl;
//            g[e].zero_one = 1;

            graph_recursive_construct(g, vertices_descriptors, v_sink_one, v_sink_zero, branch_instructions,
                    one_child_instruction_number, branch_instruction_number, 1);
        }
    } else{ // observed before
        cout << "EDGE ADDED FOR THE OBSERVED NODE #" << g[V_global[branch_instruction_number]].label <<
        " FROM BRANCH INSTRUCTION NUMBER #" << branch_instruction_number << endl;

        bool inserted = false;
        Graph::edge_descriptor e;
        Graph::vertex_descriptor parent_vertex = V_global[parent_branch_instruction_number];
        boost::tie(e, inserted) = add_edge(parent_vertex, V_global[branch_instruction_number], g);
        g[e].zero_one = zero_one;

        // updating its path count
        g[V_global[branch_instruction_number]].path_count += g[V_global[parent_branch_instruction_number]].path_count;

    }

    //TODO::DONE. correct the ordering for all add_edges so as not to mistakenly reverse the ordering of the directed graph!
    //TODO: Calculate the m_e using dynamic programming.
    //NOTE: node 1 and sink node do NOT have conflicts! Sink node key in the V_global map is 1 which corresponds to its
    // BRANCH INSTRUCTION which correctly is set to 1. Node 1, has another branch instruction number (134 for example).
    // So anywhere inside the code, V_global[1] represents the sink node which is defined at the top before graph
    // construction.


    return;
    /////////////////////////////////

}
