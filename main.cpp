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

std::map<int, Graph::vertex_descriptor> V;
//std::map<int, Graph::vertex_descriptor> V;


map<int, vector<int>> get_input(std::string filename, std::vector<std::vector<int>>& instances, std::vector<int> & labels,  int* number_of_samples, int* number_of_base_hypotheses);
void graph_recursive_construct(Graph& g, vector<Graph::vertex_descriptor>& vertices, vector<Graph::vertex_descriptor>& v_sink_one, vector<Graph::vertex_descriptor>& v_sink_zero, map<int, vector<int>> branch_instructions, int branch_instruction_number);

int main(int argc, char *argv[]){


    int number_of_samples, number_of_base_hypotheses;
    int i, j, k;
    double eps =0.0001;
    std::vector<std::vector<int>> instances;
    std::vector<int> labels;
    std::vector<int>::iterator itr;
    //double xvals[];

//    std::string filename = argv[1];
    std::string filename = "rofk-100-p.kf";



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
    graph_recursive_construct(g, V,V_sink_one_edge, V_sink_zero_edge, branch_instructions, last_instruction_number); // Fourth element is the pointer to the desired

    // when the above command is executed, graph is constructed along with its one and zero edges. Just edges to the
    // sink node remain to be added. next step is to do the topological sorting. BUT I don't think that the sorting is
    // required, since when kf format file is read from the end, the resulting recursive graph constructed, is itself
    // sorted.
    //NOTE: So we only need to add one-edge and zero-edge of the nodes in the V_sink_one_edge and V_sink_zero_edge
    // vectors. Then transform the graph to the NZDD form with a single label on each edge.

    // adding the sink node
    auto sink_node = add_vertex(g);
    // Adding the one-edges
    for (auto first = V_sink_one_edge.begin(), last = V_sink_one_edge.end(); first != last; ++first) {
        cout << "ONE-EDGE to the sink node was added. From " << (g[*first]).label << endl;

        bool inserted = false;
        Graph::edge_descriptor e;

        boost::tie(e, inserted) = add_edge(*first, sink_node, g);
        g[e].zero_one = 1;
    }
    for (auto first = V_sink_zero_edge.begin(), last = V_sink_zero_edge.end(); first != last; ++first) {
        cout << "ZERO-EDGE to the sink node was added. From " << (g[*first]).label << endl;

        bool inserted = false;
        Graph::edge_descriptor e;

        boost::tie(e, inserted) = add_edge(*first, sink_node, g);
        g[e].zero_one = 0;
    }

    //SECTION: transformation of the graph to NZDD

    // In this part we have to iterate over all edges and assign their parent's label to them.
    //NOTE: Only one-edges, get their parent's label and 0 (or empty set) will be assigned to zero-edges.

//    pair<edge_iterator, edge_iterator> ei = edges(g);
//
//    cout << "Number of edges = " << num_edges(g) << "\n" << endl;
//    cout << "Edge list:\n" << endl;
//
//    copy( ei.first, ei.second, ostream_iterator<Graph::edge_descriptor>{
//                       cout, "\n"});
//
//    cout << endl;

//    pair<edge_iterator, edge_iterator> edgePair;
//    for (edgePair = edges(g); edgePair.first != edgePair.second; ++edgePair.first)
//    {
//        std::cout << g[*edgePair.first].zero_one << std::endl;
//    }

//    std::cout << "edges(g) = ";
//    boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
//    typename GraphTraits::edge_descriptor e;
//    for (boost::tie(out_i, out_end) = out_edges(v, g);
//         out_i != out_end; ++out_i) {
//        e = *out_i;
//        Vertex src = source(e, g), targ = target(e, g);
//        std::cout << "(" << index[src] << ","
//                  << index[targ] << ") ";
//    }
//    std::cout << std::endl;
//    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
//        std::cout << "(" << (source(*ei, g)) << "," << (target(*ei, g)) << ") ";
//    std::cout << std::endl;
//
//    cout << "#####################  " << (V[10]).label << endl;

//    boost::graph_traits < Graph >::out_edge_iterator ei, ei_end;
//    for (boost::tie(ei, ei_end) = out_edges(V[100], g); ei != ei_end; ++ei) {
//        auto source = boost::source ( *ei, g );
//        auto target = boost::target ( *ei, g );
//        std::cout << "There is an edge from " << V[source].label <<  " to " << target << std::endl;
//    }
//    cout << "##########\t########## " << (*V[2]).label << endl;

    //SECTION: calculating the m_e for each edge using dynamic programming

    // instruction of the instructions vector.

    /*
      for (itr =labels.begin(); itr !=labels.end(); itr++){
      cout << *itr <<endl;
      }


    for (i=0; i<number_of_samples; i++){
      cout <<"instances[" << i <<"]={";
      for (j=0; j<instances[i].size(); j++){
        cout << instances[i][j] << " ";
      }
      cout <<"}"<< endl;
    }
    */
    /*
    for (i=0; i<m_neg; i++){
      cout <<"X_neg[" << i <<"]={";
      for (j=0; j<X_neg[i].size(); j++){
        cout << X_neg[i][j] << " ";
      }
      cout <<"}"<< endl;
    }
    */



//    try{
//        GRBEnv env = GRBEnv();
//        GRBModel model = GRBModel(env);
//
//
//        //setting
//        model.set(GRB_DoubleParam_OptimalityTol, eps);
//
//
//        //variables
//        //w
//        auto w =std::vector<GRBVar>{};
//        for (j =0; j< number_of_base_hypotheses; j++){
//            w.push_back(model.addVar(0,  GRB_INFINITY, 0,  GRB_CONTINUOUS));
//        }
//
//        //b
//        auto b = model.addVar(- GRB_INFINITY,  GRB_INFINITY, 0,  GRB_CONTINUOUS);
//
//
//        model.update();
//
//        //setting objective
//        GRBLinExpr objective = 0;
//        for (j=0; j<number_of_base_hypotheses; j++){
//            objective +=  w[j];
//        }
//        model.setObjective(objective, GRB_MINIMIZE);
//        //model.addVars(-1*GRB_INFINITY, GRB_INFINITY, c, GRB_CONTINUOUS, NULL, 10);
//        //model.addVars(-1, 1, c, GRB_CONTINUOUS, NULL, 10);
//        //GRBLinExpr expr = new GRBLinExpr();
//        //expr = x[1] + x[2];
//        //expr->addTerms(c,x,10);
//        //expr->addTerms(1.0, x[0],1); //expr.addTerms(2.0, x[1],1);
//
//
//        model.update();
//
//
//        //constraints
//        auto cnstr =std::vector<GRBConstr>{};
//        GRBLinExpr lhs = 0;
//
//        // simplex constraint
//        //for (i=0; i<number_of_samples; i++){ expr +=x[i]; }
//        //cnstr.push_back(model.addConstr(expr, GRB_EQUAL, 1.0));
//
//        //constraints for instances
//        for (i=0; i< number_of_samples; i++){
//            lhs = 0;
//            for (k=0; k<instances[i].size(); k++){ lhs += labels[i] * w[instances[i][k]]; }
//            lhs += labels[i]*b;
//            cnstr.push_back(model.addConstr(lhs, GRB_GREATER_EQUAL, 1.0));
//        }
//
//        model.update();
//
//
//
//        //optimization
//        model.optimize();
//
//        //show results
//        cout << "optval=" << model.get(GRB_DoubleAttr_ObjVal) <<endl;
//
//        cout  << "w={ ";
//        for (j=0; j<number_of_base_hypotheses; j++){
//            cout << w[j].get(GRB_DoubleAttr_X) << " ";
//        }
//        cout << " }" << endl;
//
//        cout << "b=" << b.get(GRB_DoubleAttr_X) << endl;
//
//
//        //double xvals[] = model.get(GRB_DoubleAttr_X, model.getVars());
//        //double xvals[] = model.get(GRB_DoubleAttr_X);
//        /*
//        cout  << "d={ ";
//          for (i=0; i<number_of_samples; i++){
//        cout << x[i].get(GRB_DoubleAttr_X) << " ";
//          }
//        cout << " }" << endl;
//        cout <<"gamma=" << x[number_of_samples].get(GRB_DoubleAttr_X) << endl;
//        double bias = -1 * cnstr[1].get(GRB_DoubleAttr_Pi);
//        cout << "bias="  << bias <<endl;
//        std::vector<double> alpha;
//        alpha.resize(number_of_base_hypotheses);
//        for (j=0; j<number_of_base_hypotheses; j++){
//          alpha[j]= -1* cnstr[j+2].get(GRB_DoubleAttr_Pi);
//          cout << "alpha_" <<"j=" << alpha[j] <<endl;
//        }
//        */
//
//    }catch(GRBException e) {
//        cout << "Error code = " << e.getErrorCode() << endl;
//        cout << e.getMessage() << endl;
//    } catch (...) {
//        cout << "Error during optimization" << endl;
//    }
//
//    //test(filename);
//    //test(filename, labels,H, &number_of_samples);

    cout << "Program Finished!" << endl;

    return 0;

}

map<int, vector<int>> get_input(std::string filename, std::vector<vector<int>>& instances, std::vector<int> & labels,  int* number_of_samples, int* number_of_base_hypotheses){


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

void graph_recursive_construct(Graph& g, vector<Graph::vertex_descriptor>& vertices_descriptors, vector<Graph::vertex_descriptor>& v_sink_one, vector<Graph::vertex_descriptor>& v_sink_zero, map<int, vector<int>> branch_instructions, int branch_instruction_number){

    // Termination state for recursive function:
//        if(branch_instruction_number == 0 || branch_instruction_number == 1){
//            return;
//        }
//    for (vector<vector<int>>::reverse_iterator i = branch_instructions.rbegin(); i != branch_instructions.rend(); i++) {
    cout << "Node index: " << (branch_instructions[branch_instruction_number]).at(0) << endl;
    cout << "One-edge: " << (branch_instructions[branch_instruction_number]).at(1) << endl;
    cout << "Zero-edge: " << (branch_instructions[branch_instruction_number]).at(2) << endl;

//        Graph ::vertex_descriptor temporary_descriptor;
//        temporary_descriptor = add_vertex(g);
//        vertices_descriptors.push_back(temporary_descriptor);
//        g[temporary_descriptor].label = (*i).at(0);

//        cout << "Size of vertex descriptors: " << vertices_descriptors.size() << endl;

    ///////////////////////////////// Recursive

    //NOTE: if some node is observed, then we should NOT add a vertex for it. We should just connect it to the previously
    // observed node and terminate recursion through this branch since if that node is observed previously, certainly
    // recursion has been called on its children in another recursion path. For this purpose we should keep a GLOBAL map
    // of observed nodes available to all instances of recursion calls being undertaken at each moment.
    // map<branch_instruction_number, vertex_descriptor> or
    // map<is_observed, branch_instruction_number, vertex_descriptor> could be used.

    //NOTE: To calculate the m_e, at the instant of creation of a vertex, a number reflecting the number of paths from
    // the root to this point is inherited from the parent. If at any future instance of time, an observed node is the
    // child of a another node, number of that node is added to that of this node's hence it gets updated.

    Graph ::vertex_descriptor temporary_descriptor;
    temporary_descriptor = add_vertex(g);
    vertices_descriptors.push_back(temporary_descriptor);
    // Nodes label is set to hypothesis node and later will be given to its one-edge.
    g[temporary_descriptor].label = (branch_instructions[branch_instruction_number]).at(0);
//        g[temporary_descriptor].label = (*branch_instruction_pointer).at(0);
    cout << "Node " << g[temporary_descriptor].label << " was added." << endl;

    // Defining children of this node
    Graph ::vertex_descriptor temporary_descriptor_one_child,temporary_descriptor_zero_child;

//        temporary_descriptor_one_child = add_vertex(g);
//        temporary_descriptor_zero_child = add_vertex(g);

//        vertices_descriptors.push_back(temporary_descriptor_one_child);
//        vertices_descriptors.push_back(temporary_descriptor_zero_child);

//        cout << "One child and zero child temporary descriptors added to the graph and to the descriptors vector." << endl;

    /*TODO: 1. adding the edges.
     *      2. For zero edges empty edge should be considered
     *      */

    // 1. if the one-child is zero, do nothing and leave recursion and the descendant.
    // 2. if the one-child is one, connect the ONE-edge to SINK node and finish recursion. For this reason we need a vector of descriptors that
    //      need to be connected to the sink node at the end of recursion.
    // 3. otherwise continue to its descendants and continue the recursion.
    if((branch_instructions[branch_instruction_number]).at(1) == 0 || (branch_instructions[branch_instruction_number]).at(1) == 1){
        // case 1
        if ((branch_instructions[branch_instruction_number]).at(1) == 0){
            cout << "one-child is zero and recursion stops." << endl;

            //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
            // nothing should be happened according to these nodes.
        }
            // case 2
        else {
            //TODO: Add to sink_descriptors
            v_sink_one.push_back(temporary_descriptor);
            cout << "one-child is one and recursion stops. It connects to the sink node later." << endl;
            //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
            // nothing should be happened according to these nodes.
        }
        // if goes to the following else, it means that recursion must be continued for one-edge.
    }
        // case 3
    else {
        /*TODO: Pay attention to the shifts in kf data. some start from 2!
     * because of this, we use map instead of a pure vector.
     * */

        temporary_descriptor_one_child = add_vertex(g);
        vertices_descriptors.push_back(temporary_descriptor_one_child);
        // there are non-zero and non-one children, hence we go through the recursion to consider their branch instructions
        auto one_child_instruction_number = branch_instructions[branch_instruction_number].at(1);
//            auto zero_child_instruction_number = branch_instructions[branch_instruction_number].at(2);

        // This is the one-edge.
        g[temporary_descriptor_one_child].label = branch_instructions[one_child_instruction_number].at(0);

        bool inserted = false;
        Graph::edge_descriptor e;

        boost::tie(e, inserted) = add_edge(temporary_descriptor, temporary_descriptor_one_child, g);
        std::cout<< "connecting node #" << g[temporary_descriptor].label << " to the node (ONE-EDGE) #" << g[temporary_descriptor_one_child].label << std::endl;
        g[e].zero_one = 1;

        graph_recursive_construct(g, vertices_descriptors, v_sink_one, v_sink_zero, branch_instructions, one_child_instruction_number);
    }

    // 1. if the zero-child is zero, do nothing and leave recursion and the descendant.
    // 2. if the zero-child is one, connect the ZERO-edge to SINK node and finish recursion. For this reason we need a vector of descriptors that
    //      need to be connected to the sink node at the end of recursion.
    // 3. otherwise continue to its descendants and continue the recursion.
    if((branch_instructions[branch_instruction_number]).at(2) == 0 || (branch_instructions[branch_instruction_number]).at(2) == 1){
        // case 1
        if ((branch_instructions[branch_instruction_number]).at(2) == 0){
            cout << "zero-child is zero and recursion stops." << endl;
            //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
            // nothing should be happened according to these nodes.
        }
            // case 2
        else {
            //TODO: Add to sink_descriptors
            v_sink_zero.push_back(temporary_descriptor);
            cout << "zero-child is one and recursion stops. It connects to the sink node later." << endl;
            //NOTE: If entered this if, no return should be invoked! Other if and else have to be completed! Simply
            // nothing should be happened according to these nodes.
        }
        // if goes to the following else, it means that recursion must be continued for zero-edge.
    }
        // case 3
    else {
        /*TODO: Pay attention to the shifts in kf data. some start from 2!
     * because of this, we use map instead of a pure vector.
     * */

        temporary_descriptor_zero_child = add_vertex(g);
        vertices_descriptors.push_back(temporary_descriptor_zero_child);
        // there are non-zero and non-one children, hence we go through the recursion to consider their branch instructions
        auto zero_child_instruction_number = branch_instructions[branch_instruction_number].at(2);
        //            auto zero_child_instruction_number = branch_instructions[branch_instruction_number].at(2);

        // This is the zero-edge.
        g[temporary_descriptor_zero_child].label = branch_instructions[zero_child_instruction_number].at(0);

        bool inserted = false;
        Graph::edge_descriptor e;

        boost::tie(e, inserted) = add_edge(temporary_descriptor, temporary_descriptor_zero_child, g);
        std::cout<< "connecting node #" << g[temporary_descriptor].label << " to the node (ZERO-EDGE) #" << g[temporary_descriptor_zero_child].label << std::endl;
        g[e].zero_one = 0;

        graph_recursive_construct(g, vertices_descriptors, v_sink_one, v_sink_zero, branch_instructions, zero_child_instruction_number);
    }
    return;
    /////////////////////////////////

}
