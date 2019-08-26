#include "POA.h"
#include <stdlib.h>
#include "Correct.h"
#include "Process_Read.h"
#define INIT_EDGE_SIZE 50
#define INCREASE_EDGE_SIZE 5
#define INIT_NODE_SIZE 16000




void init_Edge_alloc(Edge_alloc* list)
{
    if (list->list == NULL)
    {
        list->size = INIT_EDGE_SIZE;
        list->length = 0;
        list->list = (Edge*)malloc(sizeof(Edge)*list->size);
    }
    else
    {
        list->length = 0;
    }
    
}

void clear_Edge_alloc(Edge_alloc* list)
{
    list->length = 0;
}

void destory_Edge_alloc(Edge_alloc* list)
{
    free(list->list);
}

void append_Edge_alloc(Edge_alloc* list,  uint64_t in_node, uint64_t out_node, uint64_t weight, uint64_t length)
{
    if (list->length + 1 > list->size)
    {
        list->size = list->size + INCREASE_EDGE_SIZE;
        list->list = (Edge*)realloc(list->list, sizeof(Edge)*list->size);
    }

    list->list[list->length].in_node = in_node;
    list->list[list->length].out_node = out_node;
    list->list[list->length].weight = weight;
    list->list[list->length].length = length;
    list->list[list->length].num_insertions = 0;

    list->length++;
}




void init_Node_alloc(Node_alloc* list)
{
    list->size = INIT_NODE_SIZE;
    list->length = 0;
    list->list = (Node*)malloc(sizeof(Node)*list->size);
    list->sort.size = 0;
    list->sort.list = NULL;
    list->sort.visit = NULL;

    list->sort.iterative_buffer = NULL;
    list->sort.iterative_buffer_visit = NULL;


    long long i;
    for (i = 0; i < list->size; i++)
    {
        list->list[i].insertion_edges.list=NULL;
        list->list[i].mismatch_edges.list=NULL;
        list->list[i].deletion_edges.list=NULL;
    }    
}

void destory_Node_alloc(Node_alloc* list)
{
    uint64_t i =0;
    for (i = 0; i < list->length; i++)
    {
        destory_Edge_alloc(&list->list[i].deletion_edges);
        destory_Edge_alloc(&list->list[i].insertion_edges);
        destory_Edge_alloc(&list->list[i].mismatch_edges);
    }

    free(list->list);
    free(list->sort.list);
    free(list->sort.visit);
    free(list->sort.iterative_buffer);
    free(list->sort.iterative_buffer_visit);
    ///free(list->topo_order);
}

void clear_Node_alloc(Node_alloc* list)
{
    uint64_t i =0;
    for (i = 0; i < list->length; i++)
    {
        clear_Edge_alloc(&list->list[i].insertion_edges);
        clear_Edge_alloc(&list->list[i].mismatch_edges);
        clear_Edge_alloc(&list->list[i].deletion_edges);
    }

    list->length = 0;
}


uint64_t append_Node_alloc(Node_alloc* list, char base)
{
   
    if (list->length + 1 > list->size)
    {
        long long i = list->size;

        ///list->topo_order这里用不到，所以不用先分配空间
        ///但是还是一起分配了吧，免得麻烦
        list->size = list->size * 2;
        list->list = (Node*)realloc(list->list, sizeof(Node)*list->size);
        ///list->topo_order = (uint64_t*)realloc(list->topo_order, sizeof(uint64_t)*list->size);
        
        for (; i < list->size; i++)
        {
            list->list[i].deletion_edges.list=NULL;
            list->list[i].insertion_edges.list=NULL;
            list->list[i].mismatch_edges.list=NULL;
        }
    }
    
    list->list[list->length].ID = list->length;
    list->list[list->length].base = base;
    list->list[list->length].weight = 1;
    list->list[list->length].num_insertions = 0;
    init_Edge_alloc(&list->list[list->length].deletion_edges);
    init_Edge_alloc(&list->list[list->length].insertion_edges);
    init_Edge_alloc(&list->list[list->length].mismatch_edges);
    
    list->length++;

    return list->length - 1;
}











void init_Graph(Graph* g)
{
    init_Node_alloc(&g->g_nodes);
    g->g_n_edges = 0;
    g->g_n_nodes = 0;
    g->g_next_nodeID = 0;
    g->s_end_nodeID = 0;
    g->s_start_nodeID = 0;
    g->seq = NULL;
    g->seqID = (uint64_t)-1;
}

void destory_Graph(Graph* g)
{
    destory_Node_alloc(&g->g_nodes);
}

void clear_Graph(Graph* g)
{
    clear_Node_alloc(&g->g_nodes);

    g->g_n_edges = 0;
    g->g_n_nodes = 0;
    g->g_next_nodeID = 0;
    g->s_end_nodeID = 0;
    g->s_start_nodeID = 0;
    g->seq = NULL;
    g->seqID = (uint64_t)-1;
}







void addUnmatchedSeqToGraph(Graph* g, char* g_read_seq, long long g_read_length, long long* startID, long long* endID)
{
    long long firstID, lastID, nodeID, i;
    firstID = -1;
    lastID = -1;

    if(g_read_length == 0)
        return;

    ///start node
    nodeID  = add_Node_Graph(g, 'S');
    firstID = nodeID;
    lastID = nodeID;


    for (i = 0; i < g_read_length; i++)
    {
        nodeID = add_Node_Graph(g, g_read_seq[i]);

        ////fprintf(stderr, "nodeID: %llu\n", nodeID);
        
        if (firstID == -1)
        {
            firstID = nodeID;
        }
        if (lastID != -1)
        {
            /** 
            ///0是match边
            add_Edge_Graph(g, lastID, nodeID, 0);
            **/
            ///只有match边长度是0
            ///mismatch边长度都是1
            append_Edge_alloc(&(g->g_nodes.list[lastID].mismatch_edges), lastID, nodeID, 1, 0);
        }

        lastID = nodeID; 
    }

    *startID = firstID;
    *endID = lastID;

    g->s_start_nodeID = firstID;
    g->s_end_nodeID = lastID;
    
}








void addmatchedSeqToGraph(Graph* backbone, long long currentNodeID, char* x_string, long long x_length, 
        char* y_string, long long y_length, CIGAR* cigar, long long backbone_start, long long backbone_end)
{
    
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    int last_operation = -1;

    
    
    
    ///note that node 0 is the start node
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///这种情况代表匹配和mismatch
        if (operation == 0 || operation == 1)
        {
            
            for (i = 0; i < operationLen; i++)
            {
                //backbone->g_nodes.list[currentNodeID].weight++;
                ///前面是插入，后面有可能是误配，也有可能是匹配
                add_mismatchEdge_weight(backbone, currentNodeID, y_string[y_i], last_operation);    
                x_i++;
                y_i++;
                currentNodeID++;
            }
        }///insertion
        else if (operation == 2)
        {
            ///cigar的起始和结尾不可能是2，所以这里-1没问题
            ///if (operationLen <= CORRECT_INDEL_LENGTH)
            {
                add_insertionEdge_weight(backbone, currentNodeID, y_string + y_i, operationLen);
                backbone->g_nodes.list[currentNodeID].num_insertions++;
            }
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            ///3是y缺字符（x多字符），也就是backbone多字符
            ///这个相当于在backbone对应字符处变成了‘——’
            ///因此可以用mismatch类似的方法处理
            ///if (operationLen <= CORRECT_INDEL_LENGTH)
            {
                ///add_deletion_to_backbone(backbone, &currentNodeID, operationLen);
                ///在编辑距离中，前面是个插入，后面是个删除，这种情况是不存在的
                ///为了保险要不还给他加上吧
                ///先不加
                add_deletionEdge_weight(backbone, currentNodeID, operationLen);
            }

            
            currentNodeID += operationLen;
            x_i += operationLen;
        }

        last_operation = operation;
        
        cigar_i++;
    }
    


    /**
    ///cigar的起始和结尾不可能是2
    if (cigar->C_C[0] == 2 || cigar->C_C[cigar->length - 1] == 2)
    {
        fprintf(stderr, "error\n");
    }
    

    if (x_i != x_length)
    {
        fprintf(stderr, "x_i: %d, x_length: %d\n", x_i, x_length);
    }

    if (y_i != y_length)
    {
        fprintf(stderr, "y_i: %d, y_length: %d\n", y_i, y_length);
    }
    **/
    
}


void debug_graph(Graph* g, long long backbone_length)
{
    long long i = 0;

    if (g->s_start_nodeID != 0 || g->s_end_nodeID != backbone_length)
    {
        fprintf(stderr, "error\n");
    }
    
    
    for (i = g->s_start_nodeID; i <= g->s_end_nodeID; i++)
    {
        if(g->g_nodes.list[i].weight != 1)
        {
            fprintf(stderr, "error node weight\n");
        }
        
        if(g->g_nodes.list[i].mismatch_edges.length > 4)
        {
            fprintf(stderr, "error mismatch_edges\n");
        }

        if(g->g_nodes.list[i].mismatch_edges.length < 1 && i != g->s_end_nodeID)
        {
            fprintf(stderr, "i: %d, error mismatch_edges: %d\n", i, g->g_nodes.list[i].mismatch_edges.length);
        }
    }

    
    for (i = 0; i < g->g_nodes.length; i++)
    {
        if(g->g_nodes.list[i].ID < g->s_start_nodeID || g->g_nodes.list[i].ID > g->s_end_nodeID)
        {
            if (g->g_nodes.list[i].deletion_edges.length +
                g->g_nodes.list[i].insertion_edges.length +
                g->g_nodes.list[i].mismatch_edges.length
                != 1)
            {
                fprintf(stderr, "g->s_start_nodeID: %lld\n", 
                g->s_start_nodeID);
                fprintf(stderr, "g->s_end_nodeID: %lld\n", 
                g->s_end_nodeID);
                fprintf(stderr, "deletion_edges_length: %lld\n", 
                g->g_nodes.list[i].deletion_edges.length);
                fprintf(stderr, "insertion_edges_length: %lld, \n", 
                g->g_nodes.list[i].insertion_edges.length);
                fprintf(stderr, "g->g_nodes.list[i].insertion_edges.list[0].length: %lld, \n", 
                g->g_nodes.list[i].insertion_edges.list[0].length);
                fprintf(stderr, "g->g_nodes.list[i].insertion_edges.list[1].length: %lld, \n", 
                g->g_nodes.list[i].insertion_edges.list[1].length);

                fprintf(stderr, "mismatch_edges_length: %lld\n", 
                g->g_nodes.list[i].mismatch_edges.length);
            }
            else
            {
                ///不是0肯定是1
                if (g->g_nodes.list[i].deletion_edges.length != 0)
                {
                    long long step = g->g_nodes.list[i].deletion_edges.list[0].length;
                    long long nodeID = i;
                    for (int j = 0; j < step; j++)
                    {
                        nodeID = g->g_nodes.list[nodeID].deletion_edges.list[0].out_node;
                    }
                    
                    nodeID = g->g_nodes.list[nodeID].deletion_edges.list[0].out_node;

                    if (nodeID < g->s_start_nodeID || nodeID > g->s_end_nodeID)
                    {
                        fprintf(stderr, "error\n");
                    }
                }

                if (g->g_nodes.list[i].insertion_edges.length != 0)
                {
                    
                    long long step = g->g_nodes.list[i].insertion_edges.list[0].length;
                    long long nodeID = i;


                    for (int j = 0; j < step; j++)
                    {
                        nodeID = g->g_nodes.list[nodeID].insertion_edges.list[0].out_node;
                    }

                    nodeID = g->g_nodes.list[nodeID].insertion_edges.list[0].out_node;

                    if ((nodeID < g->s_start_nodeID || nodeID > g->s_end_nodeID))
                    {
                        fprintf(stderr, "error: step: %d\n", step);
                    }
                    
                }


                if (g->g_nodes.list[i].mismatch_edges.length != 0)
                {
                    
                    long long step = g->g_nodes.list[i].mismatch_edges.list[0].length;
                    long long nodeID = i;

                    for (int j = 0; j < step; j++)
                    {
                        nodeID = g->g_nodes.list[nodeID].mismatch_edges.list[0].out_node;
                    }

                    nodeID = g->g_nodes.list[nodeID].mismatch_edges.list[0].out_node;

                    if (nodeID < g->s_start_nodeID || nodeID > g->s_end_nodeID)
                    {
                        fprintf(stderr, "error\n");
                    }
                    
                }
                
            }
            
            
        }
    }
    


    
}

void Graph_debug(Graph* backbone, long long currentNodeID, char* x_string, long long x_length, 
        char* y_string, long long y_length, CIGAR* cigar, long long backbone_start, long long backbone_end)
{
    /**
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;

    
    

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2是x缺字符（y多字符），而3是y缺字符（x多字符）
    ///while (x_i < x_len && y_i < y_len && cigar_i < cigar->length)
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///这种情况代表匹配
        if (operation == 0)
        {
            
            for (i = 0; i < operationLen; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].base != y_string[y_i])
                {
                    fprintf(stderr, "error match\n");
                }

                backbone->g_nodes.list[currentNodeID].weight--;
            
                x_i++;
                y_i++;
                currentNodeID++;
            }
        }
        else if (operation == 1)
        {
            for (i = 0; i < operationLen; i++)
            {
                if (backbone->g_nodes.list[currentNodeID].base == y_string[y_i])
                {
                    fprintf(stderr, "error mismatch 1\n");
                }

                long long mismatchID = get_alignToNode(backbone, currentNodeID, y_string[y_i]);



                if(mismatchID == -1)
                {
                    fprintf(stderr, "error mismatch 2\n");
                }
                else
                {
                    backbone->g_nodes.list[mismatchID].weight--;
                }
                

                x_i++;
                y_i++;
                currentNodeID++;
            } 
        }
        else if (operation == 2)
        {
            long long nodeID = currentNodeID - 1;
            long long mismatchID;

            for (i = 0; i < operationLen; i++)
            {
                mismatchID = get_insertion_Node(backbone, nodeID, y_string[y_i]);

                if (mismatchID == -1)
                {
                    fprintf(stderr, "error insertion 1, i: %d\n", i);
                }
                else
                {
                    backbone->g_nodes.list[mismatchID].weight--;
                }

                nodeID = mismatchID;
                
                y_i++;
            }
            ///注意这里是x_string[x_i]而不是x_string[currentNodeID]
            mismatchID = get_insertion_Node(backbone, nodeID, x_string[x_i]);
            if (mismatchID == -1)
            {
                fprintf(stderr, "error insertion 2, i: %d, x_i: %d\n", i, x_i);
            }

            
            if (mismatchID != currentNodeID)
            {
                fprintf(stderr, "error insertion 3, i: mismatchID: %d, currentNodeID: %d\n", mismatchID, currentNodeID);
            }
            
            


        }
        else if (operation == 3)
        {
            for (i = 0; i < operationLen; i++)
            {
 
                long long mismatchID = get_alignToNode(backbone, currentNodeID, 'D');

                if(mismatchID == -1)
                {
                    fprintf(stderr, "error deletion 2\n");
                }
                else
                {
                    backbone->g_nodes.list[mismatchID].weight--;
                }

                x_i++;
                currentNodeID++;
            } 
        }
        
        cigar_i++;
    }


    if (cigar->C_C[0] == 2 || cigar->C_C[cigar->length - 1] == 2)
    {
        fprintf(stderr, "error\n");
    }
    

    if (x_i != x_length)
    {
        fprintf(stderr, "x_i: %d, x_length: %d\n", x_i, x_length);
    }

    if (y_i != y_length)
    {
        fprintf(stderr, "y_i: %d, y_length: %d\n", y_i, y_length);
    }
    **/
    
}



