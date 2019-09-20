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


//NOTE: slack variables for edges are designated by beta in the algorithm draft and base hypotheses weights are denoted
// by alpha. In this code, each edge is associated with a edge_slack_variable playing the role of beta and
// base_hypothesis_weight represents alpha.


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

    float edge_slack_variable;

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

std::map<int, Graph::vertex_descriptor> V_global_positive; // key is branch instruction number
std::map<int, Graph::vertex_descriptor> V_global_negative; // key is branch instruction number
//std::map<int, Graph::vertex_descriptor> V;

// /////////////////////// GUROBI GLOBALS

//variables

// the following variable denotes the rho in the draft.
GRBVar optimal_soft_margin; // rho

// this maps keeps base hypotheses associated weights. Key is the index of base hypothesis and the value is its
// weight.
map<int , GRBVar> base_hypotheses_weights;

// this map keeps the slack variables for all of edges. Regardless of the edge being a member of positive labeled data
// or negative labeled data.
map<Graph::edge_descriptor, GRBVar> edge_slack_variables;

// map for z variables and their corresponding vertex descriptors. Regardless of the node being a member of positive labeled data
// or negative labeled data.
map<Graph::vertex_descriptor, GRBVar> weight_distance_to_root;

// nu (coefficient of slack variables multiplication for subtraction from rho)
float nu = 0.1;

// objective of the problem used everywhere in the code
GRBLinExpr objective;

//constraints
auto constraints =std::vector<GRBConstr>{};

// GUROBI GLOBALS ///////////////////////

map<int, vector<int>> get_input(std::string filename, std::vector<std::vector<int>>& instances, std::vector<int> & labels,  int* number_of_positive_samples, int* number_of_base_hypotheses);

void graph_recursive_construct_positive(Graph& g, vector<Graph::vertex_descriptor>& vertices,
        map<int, vector<int>> branch_instructions, int branch_instruction_number, int parent_branch_instruction_number,
        int zero_one, GRBModel& grb_model);
void graph_recursive_construct_negative(Graph& g, vector<Graph::vertex_descriptor>& vertices,
        map<int, vector<int>> branch_instructions, int branch_instruction_number, int parent_branch_instruction_number,
        int zero_one, GRBModel& grb_model);

void test_model(GRBModel& grb_model);

int main(int argc, char *argv[]){


    int number_of_positive_samples, number_of_negative_samples, number_of_base_hypotheses_positive,
    number_of_base_hypotheses_negative;
    int number_of_base_hypotheses;
    int i, j, k;
    double eps =0.0001;
    std::vector<std::vector<int>> instances;
    std::vector<int> labels;
    std::vector<int>::iterator itr;
    //double xvals[];

//    std::string filename = argv[1];
//    std::string filename = "small.kf";
    std::string filename_positive_data = "rofk-100-p.kf";
    std::string filename_negative_data = "rofk-100-n.kf";



    //get input and store constraints in C_pos and C_neg
    map<int, vector<int>> branch_instructions_positive;
    map<int, vector<int>> branch_instructions_negative;
    branch_instructions_positive = get_input(filename_positive_data, instances, labels, &number_of_positive_samples, &number_of_base_hypotheses_positive);
    cout << "Size of branch_instructions for positive labeled data: " << branch_instructions_positive.size() << endl;
    cout <<"number_of_positive_samples==" << number_of_positive_samples << endl;
    cout <<"number_of_base_hypotheses_positive=" << number_of_base_hypotheses_positive << endl;
    branch_instructions_negative = get_input(filename_negative_data, instances, labels, &number_of_negative_samples, &number_of_base_hypotheses_negative);
    cout << "Size of branch_instructions for positive labeled data:" << branch_instructions_negative.size() << endl;
    cout <<"number_of_negative_samples==" << number_of_negative_samples << endl;
    cout <<"number_of_base_hypotheses_negative=" << number_of_base_hypotheses_negative << endl;

    number_of_base_hypotheses = std::max(number_of_base_hypotheses_positive, number_of_base_hypotheses_negative);

    //SECTION: constructing the graph

    Graph g;
    vector<Graph::vertex_descriptor> V;

    // FOR POSITIVE DATA
    auto iter_positive = branch_instructions_positive.end();
    iter_positive--;
    int last_instruction_number_positive = (iter_positive->first);
    cout << "###### Last Instruction number for positive data:\t" << (iter_positive->first) << endl;

    // adding sink (leaf) node for positive data
    Graph::vertex_descriptor one_sink_node_positive = add_vertex(g);
    g[one_sink_node_positive].label = 0; // so as to avoid conflict with index of base hypotheses
    g[one_sink_node_positive].parent_branch_instruction = 1; // i.e. it's the sink node and parent branch instruction for it,
    // is certainly 1!

    auto result_positive = V_global_positive.insert(std::pair<int, Graph::vertex_descriptor>(1, one_sink_node_positive));

    // FOR NEGATIVE DATA
    auto iter_negative = branch_instructions_negative.end();
    iter_negative--;
    int last_instruction_number_negative = (iter_negative->first);
    cout << "###### Last Instruction number for negative data:\t" << (iter_negative->first) << endl;

    // adding sink (leaf) node for positive data
    Graph::vertex_descriptor one_sink_node_negative = add_vertex(g);
    g[one_sink_node_negative].label = 0; // so as to avoid conflict with index of base hypotheses
    g[one_sink_node_negative].parent_branch_instruction = 1; // i.e. it's the sink node and parent branch instruction for it,
    // is certainly 1!

    auto result_negative = V_global_negative.insert(std::pair<int, Graph::vertex_descriptor>(1, one_sink_node_negative));



    /// *************************************************

    // recursively constructing the graph

    // when the above command is executed, graph is constructed along with its one and zero edges. Next step is to do
    // the topological sorting. BUT I don't think that the sorting is required, since when kf format file is read from
    // the end, the resulting recursive graph constructed, is itself sorted.
    //NOTE: So we only need to add one-edge and zero-edge of the nodes in the V_sink_one_edge and V_sink_zero_edge
    // vectors. Then transform the graph to the NZDD form with a single label on each edge.

    //SECTION: transformation of the graph to NZDD

    // In this part we have to iterate over all edges and assign their parent's label to them.
    //NOTE: Only one-edges, get their parent's label and 0 (or empty set) will be assigned to zero-edges.

    //SECTION: calculating the m_e for each edge using dynamic programming.
    // Actually the calculation is done implicitly in the graph construction and we just print it here.

    /// *************************************************

    //SECTION: LP problem set-up

    //TODO:
    // Done ---- 1. to calculate nu, we need the number of samples to multiply by some percentage. How to get the number
    //      of samples from a ZDD structure?
    // Done ---- 2. is there a difference between zero or one edges? Why were we cautious about them in the first place?
    // Done ---- 3. How to get the z_leaf? How to get z since they are minimum of a weight sum, coming from different paths.
    //      Shall I use dynamic programming?
    // Done ---- 4. I would need a map with key being source and target nodes and value being the corresponding edge descriptor.
    // Done ---- 5. In addition there's also a need for incorporating a map of which keys are base hypothesis and its
    //      keys are their weights.
    // Done ---- 5. Incorporate both positive and negative labeled data for this part.

    // note: Also slack variables should be declared and included in the constraints upon the recursive construction
    //  of the graph.


    // ------------------------------------------------------------------------------------------
    try{
        // Create new environment
        GRBEnv env = GRBEnv();

        env.set("LogFile", "mip1.log");

        // starting the environment
        env.start();

        // create an empty optimization model
        GRBModel model = GRBModel(env);

        // setting feasibility tolerance
        model.set(GRB_DoubleParam_OptimalityTol, eps);

        // variables

        // rho
        optimal_soft_margin = model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);

        GRBLinExpr b_hypotheses_constraint = 0;

        // base_hypotheses_weights
        for (int j=1; j<=number_of_base_hypotheses; ++j){
            base_hypotheses_weights[j] = model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
            b_hypotheses_constraint += base_hypotheses_weights[j];
        }

        weight_distance_to_root[V_global_positive[1]] = model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
        weight_distance_to_root[V_global_negative[1]] = model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);

        //bias
        auto bias = model.addVar(- GRB_INFINITY,  GRB_INFINITY, 0,  GRB_CONTINUOUS);

        model.update();

        //setting objective
        objective = 0;

        // adding the rho to the objective
        objective += optimal_soft_margin;

        //note: completing the objective function.

        // base hypotheses weights should sum up to one
        constraints.push_back(model.addConstr(b_hypotheses_constraint, GRB_EQUAL, 1));

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        graph_recursive_construct_positive(g, V, branch_instructions_positive, last_instruction_number_positive, 0, 0, model); // Fourth element is the pointer to the desired
        graph_recursive_construct_negative(g, V, branch_instructions_negative, last_instruction_number_negative, 0, 0, model); // Fourth element is the pointer to the desired

        model.update();

        // constraint for the leaf node (sink node), after graph construction all the elements are in the model and can be used.
        constraints.push_back(model.addConstr(weight_distance_to_root[V_global_positive[1]], GRB_GREATER_EQUAL, optimal_soft_margin));
        constraints.push_back(model.addConstr(weight_distance_to_root[V_global_negative[1]], GRB_GREATER_EQUAL, optimal_soft_margin));

        // Actually the calculation is done implicitly in the graph construction and we just print it here.
        std::pair<boost::adjacency_list<>::vertex_iterator,
        boost::adjacency_list<>::vertex_iterator> vs = boost::vertices(g);

        for (boost::adjacency_list<>::vertex_iterator v_itr = vs.first; v_itr != vs.second; ++v_itr){
            cout << "m_e for node #" << g[(*v_itr)].label << " from parent branch instruction #" <<
                 g[(*v_itr)].parent_branch_instruction << " is: " << g[(*v_itr)].path_count << endl;
        }

        for (auto ep = edges(g); ep.first != ep.second; ++ep.first)
        {
            // Get the two vertices that are joined by this edge.
            auto u = source(*ep.first,g);
            auto v = target(*ep.first,g);

            objective -= (1/nu)*((g[u].path_count)*edge_slack_variables[*ep.first]);
        }

        model.setObjective(objective, GRB_MAXIMIZE);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        model.update();


        cout << "############################### GUROBI LOG ###############################" << endl;

        //optimization
        model.optimize();

        cout << "############################### GUROBI LOG ###############################" << endl;

        //show results
        cout << "optval=" << model.get(GRB_DoubleAttr_ObjVal) <<endl;

        cout  << "base_hypotheses_weights={ ";
        for (j=1; j<=number_of_base_hypotheses; j++){
            cout << base_hypotheses_weights[j].get(GRB_DoubleAttr_X) << " ";
        }
        cout << " }" << endl;

        cout << "-----------------------------" << endl;

        cout  << "weight_distance_to_roots_positive={ ";
        for (int index=2; index <= number_of_positive_samples; ++index){
            cout << (weight_distance_to_root[V_global_positive[index]]).get(GRB_DoubleAttr_X) << " ";
        }
        cout << " }" << endl;

        cout  << "weight_distance_to_roots_negative={ ";
        for (int index=2; index <= number_of_negative_samples; ++index){
            cout << (weight_distance_to_root[V_global_negative[index]]).get(GRB_DoubleAttr_X) << " ";
        }
        cout << " }" << endl;

        cout << "-----------------------------" << endl;

        cout  << "slack_variables={ ";
        for (auto ep = edges(g); ep.first != ep.second; ++ep.first)
        {
            cout << edge_slack_variables[*ep.first].get(GRB_DoubleAttr_X) << " ";
        }

        cout << " }" << endl;

        cout << "-----------------------------" << endl;

        cout << "bias=" << bias.get(GRB_DoubleAttr_X) << endl;

    }catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during optimization" << endl;
    }

    // ------------------------------------------------------------------------------------------

    //SECTION: Done :) --------
    cout << "Program Finished!" << endl;
    cout << number_of_base_hypotheses << endl;
    cout << number_of_positive_samples << endl;
    cout << number_of_negative_samples << endl;

    return 0;

}

map<int, vector<int>> get_input(std::string filename, std::vector<vector<int>>& instances, std::vector<int> & labels,
        int* number_of_positive_samples, int* number_of_base_hypotheses){


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


    *number_of_positive_samples = index_of_sample;
    *number_of_base_hypotheses = dimensionality;

    return branch_instructions;
}

void graph_recursive_construct_positive(Graph& g, vector<Graph::vertex_descriptor>& vertices_descriptors,
        map<int, vector<int>> branch_instructions, int branch_instruction_number, int parent_branch_instruction_number,
        int zero_one, GRBModel& grb_model){

    cout << "Current branch instruction number: " << branch_instruction_number << "\tNode index: " <<
    (branch_instructions[branch_instruction_number]).at(0) << "\tZero-edge: " <<
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
    if(V_global_positive.find(branch_instruction_number) == V_global_positive.end()){ // not observed before
        cout << "Unobserved node" << endl;
        Graph ::vertex_descriptor temporary_descriptor;
        temporary_descriptor = add_vertex(g);
        vertices_descriptors.push_back(temporary_descriptor);
        g[temporary_descriptor].label = (branch_instructions[branch_instruction_number]).at(0);
        g[temporary_descriptor].parent_branch_instruction = branch_instruction_number; // if it's a root it will be zero,
        // otherwise is also correct

        V_global_positive.insert(std::pair<int, Graph ::vertex_descriptor>(branch_instruction_number, temporary_descriptor));

        // Nodes label is set to hypothesis node and later will be given to its one-edge.
        cout << "486 --- Node " << g[temporary_descriptor].label << " for branch instruction #" << branch_instruction_number <<
        " was added." << endl;

        // connecting this node to its parent and inheriting the path count from the parent and adding the corresponding
        // constraint.
        if (parent_branch_instruction_number == 0){// i.e. the node is a root itself and no edge to parent is required.
            // path count is initialized to 1.
            g[V_global_positive[branch_instruction_number]].path_count = 1;

            //note: if the following peace of code is executed before the node is created, It will become an observed node and
            // the whole program collapses!

            // for root, the weight distance to the root is zero.
            //todo: Is this correct?
            weight_distance_to_root[V_global_positive[branch_instruction_number]] = grb_model.addVar(0,  0, 0,  GRB_CONTINUOUS);
            grb_model.update();

        } else{ // i.e. the node is not a root and needs to be connected to its parent via a one or zero edge and
            // inherit the path count of its parent.
            //NOTE: Only add and inherit the path_count from the parent if they are connected by a one-edge.
            //TODO: What to do for zero-edges? For now I'll add for both zero and one edges since I don't see any
            // difference
            //NOTE: Also at this stage, constraint will be added for this pair of nodes and edge connecting them to
            // each other.

            // adding its z variable
            weight_distance_to_root[V_global_positive[branch_instruction_number]] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
            grb_model.update();

            bool inserted = false;
            Graph::edge_descriptor e;

            auto parent_descriptor = V_global_positive[parent_branch_instruction_number];
            boost::tie(e, inserted) = add_edge(parent_descriptor, temporary_descriptor, g);
            g[e].zero_one = zero_one;
            edge_slack_variables[e] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
            grb_model.update();

            g[V_global_positive[branch_instruction_number]].path_count = g[V_global_positive[parent_branch_instruction_number]].path_count;

            ///////////////// Constraint

            // constraint for this two node
            GRBLinExpr lhs = weight_distance_to_root[V_global_positive[branch_instruction_number]];
            // todo: for negative labeled data there's a slight change. (plus changes to minus)

            if(zero_one == 1){
                GRBLinExpr rhs = weight_distance_to_root[V_global_positive[parent_branch_instruction_number]] + edge_slack_variables[e] +
                        base_hypotheses_weights[(branch_instructions[parent_branch_instruction_number]).at(0)];

                constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
                cout << "------------------------Unobserved node ONE-EDGE constraint between node and its parent" << endl;

            } else{
                GRBLinExpr rhs = weight_distance_to_root[V_global_positive[parent_branch_instruction_number]] + edge_slack_variables[e];

                constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
                cout << "------------------------Unobserved node ZERO-EDGE constraint between node and its parent" << endl;

            }
            grb_model.update();
        }

        // Defining children of this node. Since it had not been observed, its children have not been observed, too!

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

                // Adding the zero-edge to the sink node
                cout << "554 --- ZERO-EDGE to the sink node was added. From branch instruction #" <<
                branch_instruction_number << " , node #" << g[V_global_positive[branch_instruction_number]].label << endl;

                bool inserted = false;
                Graph::edge_descriptor e;

                boost::tie(e, inserted) = add_edge(temporary_descriptor, V_global_positive[1], g);
                g[e].zero_one = 0;
                edge_slack_variables[e] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
                grb_model.update();

                //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
                // nothing should be happened according to these nodes.

                ///////////////// Constraint for connection to sink node with a ZERO edge

                // constraint for these two nodes
                GRBLinExpr lhs = weight_distance_to_root[V_global_positive[1]];
                // todo: for negative labeled data there's a slight change. (plus changes to minus)

                GRBLinExpr rhs = weight_distance_to_root[V_global_positive[branch_instruction_number]] + edge_slack_variables[e];

                constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
                cout << "CONSTRAINT FOR ZERO EDGE TO SINK NODE" << endl;
                grb_model.update();
            }
            // if goes to the following else, it means that recursion must be continued for zero-edge.
        }
            // case 3
        else {

            // there are non-zero and non-one children, hence we go through the recursion to consider their branch instructions
            auto zero_child_instruction_number = branch_instructions[branch_instruction_number].at(1);

            graph_recursive_construct_positive(g, vertices_descriptors, branch_instructions,
                    zero_child_instruction_number, branch_instruction_number, 0, grb_model);
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

                // Adding the zero-edge to the sink node
                cout << "616 --- ONE-EDGE to the sink node was added. From branch instruction #" <<
                     branch_instruction_number << " , node #" << g[V_global_positive[branch_instruction_number]].label << endl;

                bool inserted = false;
                Graph::edge_descriptor e;

                boost::tie(e, inserted) = add_edge(temporary_descriptor, V_global_positive[1], g);
                g[e].zero_one = 1;
                edge_slack_variables[e] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
                grb_model.update();

                //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
                // nothing should be happened according to these nodes.

                ///////////////// Constraint for connection to sink node with a ONE edge

                // constraint for this two node
                GRBLinExpr lhs = weight_distance_to_root[V_global_positive[1]];
                // todo: for negative labeled data there's a slight change. (plus changes to minus)

                GRBLinExpr rhs = weight_distance_to_root[V_global_positive[branch_instruction_number]] + edge_slack_variables[e] +
                        base_hypotheses_weights[(branch_instructions[branch_instruction_number]).at(0)];

                grb_model.update();
                cout << "CONSTRAINT FOR ONE EDGE TO SINK NODE" << endl;
                constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
                grb_model.update();

            }
            // if goes to the following else, it means that recursion must be continued for one-edge.
        }
            // case 3
        else {
            /*TODO::DONE Pay attention to the shifts in kf data. some start from 2!
         * because of this, we use map instead of a pure vector.
         * */

            // there are non-zero and non-one children, hence we go through the recursion to consider their branch instructions
            auto one_child_instruction_number = branch_instructions[branch_instruction_number].at(2);


            graph_recursive_construct_positive(g, vertices_descriptors, branch_instructions,
                    one_child_instruction_number, branch_instruction_number, 1, grb_model);
        }
    } else{ // observed before
        cout << "OBSERVED node" << endl;
        cout << "672 --- EDGE ADDED FOR THE OBSERVED NODE #" << g[V_global_positive[branch_instruction_number]].label <<
        " FROM BRANCH INSTRUCTION NUMBER #" << branch_instruction_number << endl;

        bool inserted = false;
        Graph::edge_descriptor e;
        Graph::vertex_descriptor parent_vertex = V_global_positive[parent_branch_instruction_number];
        boost::tie(e, inserted) = add_edge(parent_vertex, V_global_positive[branch_instruction_number], g);
        g[e].zero_one = zero_one;
        edge_slack_variables[e] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
        grb_model.update();

        // updating its path count
        g[V_global_positive[branch_instruction_number]].path_count += g[V_global_positive[parent_branch_instruction_number]].path_count;

        ///////////////// Constraint

        // constraint for this two node
        GRBLinExpr lhs = weight_distance_to_root[V_global_positive[branch_instruction_number]];
        // todo: for negative labeled data there's a slight change. (plus changes to minus)

        if(zero_one == 1){
            GRBLinExpr rhs = weight_distance_to_root[V_global_positive[parent_branch_instruction_number]] + edge_slack_variables[e] +
                    base_hypotheses_weights[(branch_instructions[parent_branch_instruction_number]).at(0)];

            constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
            cout << "***OBSERVED NODE AND ONE-EDGE***" << endl;
        } else{
            GRBLinExpr rhs = weight_distance_to_root[V_global_positive[parent_branch_instruction_number]] + edge_slack_variables[e];

            constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
            cout << "***OBSERVED NODE AND ZERO-EDGE***" << endl;
        }
        grb_model.update();

    }

    //TODO::DONE. correct the ordering for all add_edges so as not to mistakenly reverse the ordering of the directed graph!
    //TODO: Calculate the m_e using dynamic programming.
    //NOTE: node 1 and sink node do NOT have conflicts! Sink node key in the V_global_positive map is 1 which corresponds to its
    // BRANCH INSTRUCTION which correctly is set to 1. Node 1, has another branch instruction number (134 for example).
    // So anywhere inside the code, V_global_positive[1] represents the sink node which is defined at the top before graph
    // construction.

    cout << "One recursive call returned" << endl;

    return;
    /////////////////////////////////

}

void graph_recursive_construct_negative(Graph& g, vector<Graph::vertex_descriptor>& vertices_descriptors,
                                        map<int, vector<int>> branch_instructions, int branch_instruction_number, int parent_branch_instruction_number,
                                        int zero_one, GRBModel& grb_model){

    cout << "Current branch instruction number: " << branch_instruction_number << "\tNode index: " <<
         (branch_instructions[branch_instruction_number]).at(0) << "\tZero-edge: " <<
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
    if(V_global_negative.find(branch_instruction_number) == V_global_negative.end()){ // not observed before
        cout << "Unobserved node" << endl;
        Graph ::vertex_descriptor temporary_descriptor;
        temporary_descriptor = add_vertex(g);
        vertices_descriptors.push_back(temporary_descriptor);
        g[temporary_descriptor].label = (branch_instructions[branch_instruction_number]).at(0);
        g[temporary_descriptor].parent_branch_instruction = branch_instruction_number; // if it's a root it will be zero,
        // otherwise is also correct

        V_global_negative.insert(std::pair<int, Graph ::vertex_descriptor>(branch_instruction_number, temporary_descriptor));

        // Nodes label is set to hypothesis node and later will be given to its one-edge.
        cout << "486 --- Node " << g[temporary_descriptor].label << " for branch instruction #" << branch_instruction_number <<
             " was added." << endl;

        // connecting this node to its parent and inheriting the path count from the parent and adding the corresponding
        // constraint.
        if (parent_branch_instruction_number == 0){// i.e. the node is a root itself and no edge to parent is required.
            // path count is initialized to 1.
            g[V_global_negative[branch_instruction_number]].path_count = 1;

            //note: if the following peace of code is executed before the node is created, It will become an observed node and
            // the whole program collapses!

            // for root, the weight distance to the root is zero.
            //todo: Is this correct?
            weight_distance_to_root[V_global_negative[branch_instruction_number]] = grb_model.addVar(0,  0, 0,  GRB_CONTINUOUS);
            grb_model.update();

        } else{ // i.e. the node is not a root and needs to be connected to its parent via a one or zero edge and
            // inherit the path count of its parent.
            //NOTE: Only add and inherit the path_count from the parent if they are connected by a one-edge.
            //TODO: What to do for zero-edges? For now I'll add for both zero and one edges since I don't see any
            // difference
            //NOTE: Also at this stage, constraint will be added for this pair of nodes and edge connecting them to
            // each other.

            // adding its z variable
            weight_distance_to_root[V_global_negative[branch_instruction_number]] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
            grb_model.update();

            bool inserted = false;
            Graph::edge_descriptor e;

            auto parent_descriptor = V_global_negative[parent_branch_instruction_number];
            boost::tie(e, inserted) = add_edge(parent_descriptor, temporary_descriptor, g);
            g[e].zero_one = zero_one;
            edge_slack_variables[e] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
            grb_model.update();

            g[V_global_negative[branch_instruction_number]].path_count = g[V_global_negative[parent_branch_instruction_number]].path_count;

            ///////////////// Constraint

            // constraint for this two node
            GRBLinExpr lhs = weight_distance_to_root[V_global_negative[branch_instruction_number]];

            if(zero_one == 1){
                GRBLinExpr rhs = weight_distance_to_root[V_global_negative[parent_branch_instruction_number]] + edge_slack_variables[e] -
                                 base_hypotheses_weights[(branch_instructions[parent_branch_instruction_number]).at(0)];

                constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
                cout << "------------------------Unobserved node ONE-EDGE constraint between node and its parent" << endl;

            } else{
                GRBLinExpr rhs = weight_distance_to_root[V_global_negative[parent_branch_instruction_number]] + edge_slack_variables[e];

                constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
                cout << "------------------------Unobserved node ZERO-EDGE constraint between node and its parent" << endl;

            }
            grb_model.update();
        }

        // Defining children of this node. Since it had not been observed, its children have not been observed, too!

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
                // Connecting the zero-edge to the sink node
                cout << "554 --- ZERO-EDGE to the sink node was added. From branch instruction #" <<
                     branch_instruction_number << " , node #" << g[V_global_negative[branch_instruction_number]].label << endl;

                bool inserted = false;
                Graph::edge_descriptor e;

                boost::tie(e, inserted) = add_edge(temporary_descriptor, V_global_negative[1], g);
                g[e].zero_one = 0;
                edge_slack_variables[e] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
                grb_model.update();

                //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
                // nothing should be happened according to these nodes.

                ///////////////// Constraint for connection to sink node with a ZERO edge

                // constraint for these two nodes
                GRBLinExpr lhs = weight_distance_to_root[V_global_negative[1]];

                GRBLinExpr rhs = weight_distance_to_root[V_global_negative[branch_instruction_number]] + edge_slack_variables[e];

                constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
                cout << "CONSTRAINT FOR ZERO EDGE TO SINK NODE" << endl;
                grb_model.update();
            }
            // if goes to the following else, it means that recursion must be continued for zero-edge.
        }
            // case 3
        else {

            // there are non-zero and non-one children, hence we go through the recursion to consider their branch instructions
            auto zero_child_instruction_number = branch_instructions[branch_instruction_number].at(1);

            graph_recursive_construct_negative(g, vertices_descriptors, branch_instructions,
                                      zero_child_instruction_number, branch_instruction_number, 0, grb_model);
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
                // Connecting the zero-edge to the sink node
                cout << "616 --- ONE-EDGE to the sink node was added. From branch instruction #" <<
                     branch_instruction_number << " , node #" << g[V_global_negative[branch_instruction_number]].label << endl;

                bool inserted = false;
                Graph::edge_descriptor e;

                boost::tie(e, inserted) = add_edge(temporary_descriptor, V_global_negative[1], g);
                g[e].zero_one = 1;
                edge_slack_variables[e] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
                grb_model.update();

                //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
                // nothing should be happened according to these nodes.

                ///////////////// Constraint for connection to sink node with a ONE edge

                // constraint for this two node
                GRBLinExpr lhs = weight_distance_to_root[V_global_negative[1]];

                GRBLinExpr rhs = weight_distance_to_root[V_global_positive[branch_instruction_number]] + edge_slack_variables[e] -
                                 base_hypotheses_weights[(branch_instructions[branch_instruction_number]).at(0)];

                grb_model.update();
                cout << "CONSTRAINT FOR ONE EDGE TO SINK NODE" << endl;
                constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
                grb_model.update();

            }
            // if goes to the following else, it means that recursion must be continued for one-edge.
        }
            // case 3
        else {
            /*TODO::DONE Pay attention to the shifts in kf data. some start from 2!
         * because of this, we use map instead of a pure vector.
         * */

            // there are non-zero and non-one children, hence we go through the recursion to consider their branch instructions
            auto one_child_instruction_number = branch_instructions[branch_instruction_number].at(2);


            graph_recursive_construct_negative(g, vertices_descriptors, branch_instructions,
                                      one_child_instruction_number, branch_instruction_number, 1, grb_model);
        }
    } else{ // observed before
        cout << "OBSERVED node" << endl;
        cout << "672 --- EDGE ADDED FOR THE OBSERVED NODE #" << g[V_global_negative[branch_instruction_number]].label <<
             " FROM BRANCH INSTRUCTION NUMBER #" << branch_instruction_number << endl;

        bool inserted = false;
        Graph::edge_descriptor e;
        Graph::vertex_descriptor parent_vertex = V_global_negative[parent_branch_instruction_number];
        boost::tie(e, inserted) = add_edge(parent_vertex, V_global_negative[branch_instruction_number], g);
        g[e].zero_one = zero_one;
        edge_slack_variables[e] = grb_model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
        grb_model.update();

        // updating its path count
        g[V_global_negative[branch_instruction_number]].path_count += g[V_global_negative[parent_branch_instruction_number]].path_count;

        ///////////////// Constraint

        // constraint for this two node
        GRBLinExpr lhs = weight_distance_to_root[V_global_negative[branch_instruction_number]];

        if(zero_one == 1){
            GRBLinExpr rhs = weight_distance_to_root[V_global_positive[parent_branch_instruction_number]] + edge_slack_variables[e] -
                             base_hypotheses_weights[(branch_instructions[parent_branch_instruction_number]).at(0)];

            constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
            cout << "***OBSERVED NODE AND ONE-EDGE***" << endl;
        } else{
            GRBLinExpr rhs = weight_distance_to_root[V_global_positive[parent_branch_instruction_number]] + edge_slack_variables[e];

            constraints.push_back(grb_model.addConstr(lhs, GRB_LESS_EQUAL, rhs));
            cout << "***OBSERVED NODE AND ZERO-EDGE***" << endl;
        }
        grb_model.update();

    }

    //TODO::DONE. correct the ordering for all add_edges so as not to mistakenly reverse the ordering of the directed graph!
    //TODO: Calculate the m_e using dynamic programming.
    //NOTE: node 1 and sink node do NOT have conflicts! Sink node key in the V_global_positive map is 1 which corresponds to its
    // BRANCH INSTRUCTION which correctly is set to 1. Node 1, has another branch instruction number (134 for example).
    // So anywhere inside the code, V_global_positive[1] represents the sink node which is defined at the top before graph
    // construction.

    cout << "One recursive call returned" << endl;

    return;
    /////////////////////////////////

}

void test_model(GRBModel& model){
    GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
    GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
    GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

    // Add constraint: x + 2 y + 3 z <= 4
    model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

    // Set objective: maximize x + y + 2 z
    model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

    // Add constraint: x + y >= 1
    model.addConstr(x + y >= 1, "c1");


    model.optimize();

    cout << x.get(GRB_StringAttr_VarName) << " "
         << x.get(GRB_DoubleAttr_X) << endl;
    cout << y.get(GRB_StringAttr_VarName) << " "
         << y.get(GRB_DoubleAttr_X) << endl;
    cout << z.get(GRB_StringAttr_VarName) << " "
         << z.get(GRB_DoubleAttr_X) << endl;

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
    return;
}