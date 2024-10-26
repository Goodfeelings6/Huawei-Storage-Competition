#include "LKH.h"

/*
 * The Distance_SPECIAL function may be used to specify a user defined
 * distance fuction. The function is used when the EDGE_WEIGHT_TYPE is
 * SPECIAL. 
 * 
 * Example:
 *  
 *      int Distance_SPECIAL(Node * Na, Node * Nb) 
 *      {
 *           double dx = Na->X - Nb->X;
 *           double dy = Na->Y - Nb->Y;
 *           return (int) (1000 * sqrt(dx * dx + dy * dy) + 0.5);
 *      }           
 */
#define INF INT_MAX/200
int Distance_SPECIAL(Node * Na, Node * Nb)
{
    int len = DimensionSaved;
    const InputParam *input = OriginInput;
    // 由于是ATSP,内部的结点是做了转换的,所以求距离也得做相应转换
    if ((Na->Id <= len) == (Nb->Id <= len))
        return INF;
    else if (abs(Na->Id - Nb->Id) == len)
        return 0;
    else{
        int i, j;
        if (Na->Id <= len){
            i = Na->Id - 1;
            j = Nb->Id - 1 - len;
        }else{
            i = Nb->Id - 1;
            j = Na->Id - 1 - len;
        }
        // 以下正常对应邻接矩阵位置求距离
        if(i == j)
            return 0;
        else if(i==len-1){ // 磁头节点
            if(j==len-2){ // 磁头节点到虚拟节点的代价为极大值
                return INF;
            }
            else{ // 磁头节点到其他节点的代价为对应寻址时间
                HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
                HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
                return SeekTimeCalculate(&start, &end);  
            }
        }
        else if(i==len-2){ // 虚拟节点
            if(j==len-1){ // 虚拟节点到磁头节点的代价为0
                return 0;
            }
            else{ // 虚拟节点到其他节点的代价为极大值
                return INF;
            }
        }
        else{ // 其他节点
            if(j==len-1){// 其他节点到磁头节点的代价为极大值
                return INF;
            }
            else if(j==len-2){// 其他节点到虚拟节点的代价为0
                return 0;
            }
            else{
                HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
                HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
                return SeekTimeCalculate(&start, &end);
            }
        }
    }
}
