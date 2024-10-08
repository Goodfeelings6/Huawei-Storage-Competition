#include "LKH.h"
#include "Genetic.h"
#include "Heap.h"
#include "LKHInterface.h"

/*
 * LKH内核 接口，欲调用LHK，按以下三个步骤走：
 * (1)设定 LKH内核 的相关参数值：函数 ReadParameters()。一般情况下无需手动设置，采用默认值即可;
 * (2)设定问题类型、规模等，并进行样例数据的输入：函数 ReadProblem()。此处是ATSP问题，邻接矩阵形式输入;
 * (3)调用 LKH内核 求解：LKHmain().
 *
 * 求解完成后：
 * (1)提取求解结果以返回：函数 OutputTourResult() ;
 * (2)重置 LKH内核 状态：函数 ReSetLKH()。这一步确保当前调用不会对后续调用产生影响.
 */

static void CheckSpecificationPart(void);
static char *Copy(char *S);
static void CreateNodes(void);
static int FixEdge(Node * Na, Node * Nb);
static void Read_FIXED_EDGES_SECTION(LKHInput* lkhInput);
static void Read_TOUR_SECTION(LKHInput* lkhInput);
static void Read_EDGE_WEIGHT_SECTION(LKHInput* lkhInput);
static int TwoDWeightType(void);
static int ThreeDWeightType(void);

void loadDefaultParam(LKHParameters* p){
    p->ProblemFileName = 0;
    p->PiFileName = 0;
    p->InputTourFileName = 0;
    p->OutputTourFileName = 0;
    p->TourFileName = 0;
    p->CandidateFiles = 0;
    p->MergeTourFiles = 0;

    p->AscentCandidates = 50;
    p->BackboneTrials = 0;
    p->Backtracking = 0;
    p->BWTSP_B = 0;
    p->BWTSP_Q = 0;
    p->BWTSP_L = INT_MAX;
    p->CandidateSetSymmetric = 0;
    p->CandidateSetType = ALPHA;
    p->Crossover = ERXT;
    p->DelaunayPartitioning = 0;
    p->DelaunayPure = 0;
    p->DemandDimension = 1;
    p->DistanceLimit = DBL_MAX;
    p->Excess = -1;
    p->ExternalSalesmen = 0;
    p->ExtraCandidates = 0;
    p->ExtraCandidateSetSymmetric = 0;
    p->ExtraCandidateSetType = QUADRANT;
    p->Gain23Used = 1;
    p->GainCriterionUsed = 1;
    p->GridSize = 1000000.0;
    p->InitialPeriod = -1;
    p->InitialStepSize = 0;
    p->InitialTourAlgorithm = WALK;
    p->InitialTourFraction = 1.0;
    p->KarpPartitioning = 0;
    p->KCenterPartitioning = 0;
    p->KMeansPartitioning = 0;
    p->Kicks = 1;
    p->KickType = 0;
    p->MaxBreadth = INT_MAX;
    p->MaxCandidates = 5;
    p->MaxPopulationSize = 0;
    p->MaxSwaps = -1;
    p->MaxTrials = -1;
    p->MoorePartitioning = 0;
    p->MoveType = 5;
    p->MoveTypeSpecial = 0;
    p->MTSPDepot = 1;
    p->MTSPMinSize = 1;
    p->MTSPMaxSize = -1;
    p->MTSPObjective = -1;
    p->NonsequentialMoveType = -1;
    p->Optimum = MINUS_INFINITY;
    p->PatchingA = 1;
    p->PatchingC = 0;
    p->PatchingAExtended = 0;
    p->PatchingARestricted = 0;
    p->PatchingCExtended = 0;
    p->PatchingCRestricted = 0;
    p->Precision = 100;
    p->Probability = 100;
    p->POPMUSIC_InitialTour = 0;
    p->POPMUSIC_MaxNeighbors = 5;
    p->POPMUSIC_SampleSize = 10;
    p->POPMUSIC_Solutions = 50;
    p->POPMUSIC_Trials = 1;
    p->Recombination = IPT;
    p->RestrictedSearch = 1;
    p->RohePartitioning = 0;
    p->Runs = 0;
    p->Salesmen = 1;
    p->Scale = -1;
    p->Seed = 1;
    p->SierpinskiPartitioning = 0;
    p->StopAtOptimum = 1;
    p->Subgradient = 1;
    p->SubproblemBorders = 0;
    p->SubproblemsCompressed = 0;
    p->SubproblemSize = 0;
    p->SubsequentMoveType = 0;
    p->SubsequentMoveTypeSpecial = 0;
    p->SubsequentPatching = 1;
    p->TimeLimit = DBL_MAX;
    p->TotalTimeLimit = DBL_MAX;
    p->TraceLevel = 1;
    p->TSPTW_Makespan = 0;
    p->MaxMatrixDimension = 20006;
}

void ReadParameters(LKHInput* lkhInput)
{
    LKHParameters* lkhParam = lkhInput->lkhParameters;
    // 从结构体获取参数
    ProblemFileName = lkhParam->ProblemFileName;
    PiFileName = lkhParam->PiFileName;
    InputTourFileName = lkhParam->InputTourFileName;
    OutputTourFileName = lkhParam->OutputTourFileName;
    TourFileName = lkhParam->TourFileName;
    CandidateFiles = lkhParam->CandidateFiles;
    MergeTourFiles = lkhParam->MergeTourFiles;
    AscentCandidates = lkhParam->AscentCandidates;
    BackboneTrials = lkhParam->BackboneTrials;
    Backtracking = lkhParam->Backtracking;
    BWTSP_B = lkhParam->BWTSP_B;
    BWTSP_Q = lkhParam->BWTSP_Q;
    BWTSP_L = lkhParam->BWTSP_L;
    CandidateSetSymmetric = lkhParam->CandidateSetSymmetric;
    CandidateSetType = lkhParam->CandidateSetType;
    Crossover = lkhParam->Crossover;
    DelaunayPartitioning = lkhParam->DelaunayPartitioning;
    DelaunayPure = lkhParam->DelaunayPure;
    DemandDimension = lkhParam->DemandDimension;
    DistanceLimit = lkhParam->DistanceLimit;
    Excess = lkhParam->Excess;
    ExternalSalesmen = lkhParam->ExternalSalesmen;
    ExtraCandidates = lkhParam->ExtraCandidates;
    ExtraCandidateSetSymmetric = lkhParam->ExtraCandidateSetSymmetric;
    ExtraCandidateSetType = lkhParam->ExtraCandidateSetType;
    Gain23Used = lkhParam->Gain23Used;
    GainCriterionUsed = lkhParam->GainCriterionUsed;
    GridSize = lkhParam->GridSize;
    InitialPeriod = lkhParam->InitialPeriod;
    InitialStepSize = lkhParam->InitialStepSize;
    InitialTourAlgorithm = lkhParam->InitialTourAlgorithm;
    InitialTourFraction = lkhParam->InitialTourFraction;
    KarpPartitioning = lkhParam->KarpPartitioning;
    KCenterPartitioning = lkhParam->KCenterPartitioning;
    KMeansPartitioning = lkhParam->KMeansPartitioning;
    Kicks = lkhParam->Kicks;
    KickType = lkhParam->KickType;
    MaxBreadth = lkhParam->MaxBreadth;
    MaxCandidates = lkhParam->MaxCandidates;
    MaxPopulationSize = lkhParam->MaxPopulationSize;
    MaxSwaps = lkhParam->MaxSwaps;
    MaxTrials = lkhParam->MaxTrials;
    MoorePartitioning = lkhParam->MoorePartitioning;
    MoveType = lkhParam->MoveType;
    MoveTypeSpecial = lkhParam->MoveTypeSpecial;
    MTSPDepot = lkhParam->MTSPDepot;
    MTSPMinSize = lkhParam->MTSPMinSize;
    MTSPMaxSize = lkhParam->MTSPMaxSize;
    MTSPObjective = lkhParam->MTSPObjective;
    NonsequentialMoveType = lkhParam->NonsequentialMoveType;
    Optimum  = lkhParam->Optimum ;
    PatchingA = lkhParam->PatchingA;
    PatchingC = lkhParam->PatchingC;
    PatchingAExtended = lkhParam->PatchingAExtended;
    PatchingARestricted = lkhParam->PatchingARestricted;
    PatchingCExtended = lkhParam->PatchingCExtended;
    PatchingCRestricted = lkhParam->PatchingCRestricted;
    Precision = lkhParam->Precision;
    Probability = lkhParam->Probability;
    POPMUSIC_InitialTour = lkhParam->POPMUSIC_InitialTour;
    POPMUSIC_MaxNeighbors = lkhParam->POPMUSIC_MaxNeighbors;
    POPMUSIC_SampleSize = lkhParam->POPMUSIC_SampleSize;
    POPMUSIC_Solutions = lkhParam->POPMUSIC_Solutions;
    POPMUSIC_Trials = lkhParam->POPMUSIC_Trials;
    Recombination = lkhParam->Recombination;
    RestrictedSearch = lkhParam->RestrictedSearch;
    RohePartitioning = lkhParam->RohePartitioning;
    Runs = lkhParam->Runs;
    Salesmen = lkhParam->Salesmen;
    Scale = lkhParam->Scale;
    Seed = lkhParam->Seed;
    SierpinskiPartitioning = lkhParam->SierpinskiPartitioning;
    StopAtOptimum = lkhParam->StopAtOptimum;
    Subgradient = lkhParam->Subgradient;
    SubproblemBorders = lkhParam->SubproblemBorders;
    SubproblemsCompressed = lkhParam->SubproblemsCompressed;
    SubproblemSize = lkhParam->SubproblemSize;
    SubsequentMoveType = lkhParam->SubsequentMoveType;
    SubsequentMoveTypeSpecial = lkhParam->SubsequentMoveTypeSpecial;
    SubsequentPatching = lkhParam->SubsequentPatching;
    TimeLimit = lkhParam->TimeLimit;
    TotalTimeLimit = lkhParam->TotalTimeLimit;
    TraceLevel = lkhParam->TraceLevel;
    TSPTW_Makespan = lkhParam->TSPTW_Makespan;
    MaxMatrixDimension = lkhParam->MaxMatrixDimension;

    /* other */
    MergeWithTour =
        Recombination == GPX2 ? MergeWithTourGPX2 :
        Recombination == CLARIST ? MergeWithTourCLARIST :
                                   MergeWithTourIPT;
}

void ReadProblem(LKHInput* lkhInput)
{   
   FreeStructures();
    FirstNode = 0;
    C = 0;
    c = 0;

    /* (1) The specification part 定义问题类型、规模和数据输入格式 */
    Name = Copy("IOSchedule");
    Type = Copy("ATSP");
    ProblemType = ATSP;
    Asymmetric = 1;
    WeightType = EXPLICIT;
    WeightFormat = FULL_MATRIX;
    EdgeWeightType = Copy("EXPLICIT");
    EdgeWeightFormat = Copy("FULL_MATRIX");
    EdgeDataFormat = NULL;
    NodeCoordType = Copy("NO_COORDS");
    DisplayDataType = 0;
    Distance = Distance_EXPLICIT;;

    Dimension = lkhInput->matDimension;
    DimensionSaved = Dim = Dimension;
    DistanceLimit = DBL_MAX;
    Salesmen = 1;

    /* (2) The data part 根据上方定义，进行实际数据输入 */
    /* NODE_COORD_SECTION */
    // NULL
    /* EDGE_DATA_SECTION */
    // NULL
    /* FIXED_EDGES_SECTION */
    Read_FIXED_EDGES_SECTION(lkhInput);
    /* DISPLAY_DATA_SECTION */
    // NULL
    /* TOUR_SECTION */
    Read_TOUR_SECTION(lkhInput);
    /* EDGE_WEIGHT_SECTION */
    Read_EDGE_WEIGHT_SECTION(lkhInput);


    Swaps = 0;

    /* Adjust parameters 调整参数。比如某个参数的默认值依赖于问题类型或输入数据，
       在此前(ReadParameter 函数中)无法确定值是多少，仅是设定为一个默认值标记，
       所以在此处(ReadProblem 函数最后)其默认值才被实际的设定 */
    int i, j, K;
    if (Seed == 0)
        Seed = (unsigned)time(0);
    if (InitialStepSize == 0)
        InitialStepSize = 1;
    if (MaxSwaps < 0)
        MaxSwaps = Dimension;
    if (KickType > Dimension / 2)
        KickType = Dimension / 2;
    if (Runs == 0)
        Runs = 10;
    if (MaxCandidates > Dimension - 1)
        MaxCandidates = Dimension - 1;
    if (ExtraCandidates > Dimension - 1)
        ExtraCandidates = Dimension - 1;
    if (SubproblemSize >= Dimension)
        SubproblemSize = Dimension;
    else if (SubproblemSize == 0) {
        if (AscentCandidates > Dimension - 1)
            AscentCandidates = Dimension - 1;
        if (InitialPeriod < 0) {
            InitialPeriod = Dimension / 2;
            if (InitialPeriod < 100)
                InitialPeriod = 100;
        }
        if (Excess < 0)
            Excess = 1.0 / DimensionSaved * Salesmen;
        if (MaxTrials == -1)
            MaxTrials = Dimension;
        HeapMake(Dimension);
    }
    if (POPMUSIC_MaxNeighbors > Dimension - 1)
        POPMUSIC_MaxNeighbors = Dimension - 1;
    if (POPMUSIC_SampleSize > Dimension)
        POPMUSIC_SampleSize = Dimension;
    Depot = &NodeSet[MTSPDepot];
    if (ProblemType == CVRP) {
        Node *N;
        int MinSalesmen;
        if (Capacity <= 0)
            eprintf("CAPACITY not specified");
        TotalDemand = 0;
        N = FirstNode;
        do
            TotalDemand += N->Demand;
        while ((N = N->Suc) != FirstNode);
        MinSalesmen =
            TotalDemand / Capacity + (TotalDemand % Capacity != 0);
        if (Salesmen == 1) {
            Salesmen = MinSalesmen;
            if (Salesmen > Dimension)
                eprintf("CVRP: SALESMEN larger than DIMENSION");
        } else if (Salesmen < MinSalesmen)
            eprintf("CVRP: SALESMEN too small to meet demand");
        assert(Salesmen >= 1 && Salesmen <= Dimension);
        if (Salesmen == 1)
            ProblemType = TSP;
        Penalty = Penalty_CVRP;
    } else if (ProblemType == SOP || ProblemType == M1_PDTSP) {
        Constraint *Con;
        Node *Ni, *Nj;
        int n, k;
        OldDistance = Distance;
        Distance = Distance_SOP;
        if (ProblemType == M1_PDTSP) {
            for (i = 2; i < Dim; i++) {
                Ni = &NodeSet[i];
                for (k = n = 0; k < DemandDimension; k++) {
                    n = Ni->M_Demand[k];
                    if (n >= 0)
                        continue;
                    for (j = 2; j < Dim; j++) {
                        if (j == i)
                            continue;
                        Nj = &NodeSet[j];
                        if (Nj->M_Demand[k] == -n) {
                            Ni->C[j] = -1;
                            break;
                        }
                    }
                }
            }
        }
        for (j = 2; j < Dim; j++) {
            Nj = &NodeSet[j];
            for (i = 2; i < Dim; i++) {
                if (i != j && Nj->C[i] == -1) {
                    Ni = &NodeSet[i];
                    assert(Con =
                        (Constraint *)malloc(sizeof(Constraint)));
                    Con->t1 = Ni;
                    Con->t2 = Nj;
                    Con->Suc = FirstConstraint;
                    FirstConstraint = Con;
                    Con->Next = Ni->FirstConstraint;
                    Ni->FirstConstraint = Con;
                }
            }
        }
        Salesmen = 1;
        Penalty = ProblemType == SOP ? Penalty_SOP : Penalty_M1_PDTSP;
    }
    if (ProblemType == TSPTW) {
        Salesmen = 1;
        Penalty = Penalty_TSPTW;
    } else
        TSPTW_Makespan = 0;
    if (Salesmen > 1) {
        if (Salesmen > Dim)
            eprintf("Too many salesmen/vehicles (> DIMENSION)");
        MTSP2TSP();
    }
    if (ExternalSalesmen > Salesmen)
        ExternalSalesmen = Salesmen;
    if (ProblemType == ACVRP)
        Penalty = Penalty_ACVRP;
    else if (ProblemType == CCVRP)
        Penalty = Penalty_CCVRP;
    else if (ProblemType == CTSP)
        Penalty = Penalty_CTSP;
    else if (ProblemType == CVRPTW)
        Penalty = Penalty_CVRPTW;
    else if (ProblemType == MLP)
        Penalty = Penalty_MLP;
    else if (ProblemType == OVRP)
        Penalty = Penalty_OVRP;
    else if (ProblemType == PDTSP)
        Penalty = Penalty_PDTSP;
    else if (ProblemType == PDTSPF)
        Penalty = Penalty_PDTSPF;
    else if (ProblemType == PDTSPL)
        Penalty = Penalty_PDTSPL;
    else if (ProblemType == PDPTW)
        Penalty = Penalty_PDPTW;
    else if (ProblemType == ONE_PDTSP)
        Penalty = Penalty_1_PDTSP;
    else if (ProblemType == M_PDTSP)
        Penalty = Penalty_M_PDTSP;
    else if (ProblemType == M1_PDTSP)
        Penalty = Penalty_M1_PDTSP;
    else if (ProblemType == RCTVRP || ProblemType == RCTVRPTW)
        Penalty = Penalty_RCTVRP;
    else if (ProblemType == TRP)
        Penalty = Penalty_TRP;
    else if (ProblemType == TSPDL)
        Penalty = Penalty_TSPDL;
    else if (ProblemType == TSPPD)
        Penalty = Penalty_TSPPD;
    if (ProblemType == VRPB)
        Penalty = Penalty_VRPB;
    else if (ProblemType == VRPBTW)
        Penalty = Penalty_VRPBTW;
    else if (ProblemType == VRPPD)
        Penalty = Penalty_VRPPD;
    if (BWTSP_B > 0) {
        if (Penalty)
            eprintf("BWTSP not compatible with problem type %s\n", Type);
        ProblemType = BWTSP;
        free(Type);
        Type = Copy("BWTSP");
        Penalty = Penalty_BWTSP;
        if (BWTSP_L != INT_MAX)
            BWTSP_L *= Scale;
    }
    if (Penalty && (SubproblemSize > 0 || SubproblemTourFile))
        eprintf("Partitioning not implemented for constrained problems");
    Depot->DepotId = 1;
    for (i = Dim + 1; i <= DimensionSaved; i++)
        NodeSet[i].DepotId = i - Dim + 1;
    if (Dimension != DimensionSaved) {
        NodeSet[Depot->Id + DimensionSaved].DepotId = 1;
        for (i = Dim + 1; i <= DimensionSaved; i++)
            NodeSet[i + DimensionSaved].DepotId = i - Dim + 1;
    }
    if (ServiceTime != 0) {
        for (i = 1; i <= Dim; i++)
            NodeSet[i].ServiceTime = ServiceTime;
        Depot->ServiceTime = 0;
    }
    if (CostMatrix == 0 && Dimension <= MaxMatrixDimension &&
        Distance != 0 && Distance != Distance_1 && Distance != Distance_LARGE
        && Distance != Distance_LARGE && Distance != Distance_ATSP
        && Distance != Distance_MTSP && Distance != Distance_SPECIAL) {
        Node *Ni, *Nj;
        assert(CostMatrix =
            (int *)calloc((size_t)Dim * (Dim - 1) / 2, sizeof(int)));
        Ni = FirstNode->Suc;
        do {
            Ni->C =
                &CostMatrix[(size_t)(Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
            if (ProblemType != HPP || Ni->Id <= Dim)
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = Fixed(Ni, Nj) ? 0 : Distance(Ni, Nj);
            else
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = 0;
        } while ((Ni = Ni->Suc) != FirstNode);
        c = 0;
        WeightType = EXPLICIT;
    }
    if (ProblemType == TSPTW ||
        ProblemType == CVRPTW || ProblemType == VRPBTW ||
        ProblemType == PDPTW || ProblemType == RCTVRPTW) {
        M = INT_MAX / 2 / Precision;
        for (i = 1; i <= Dim; i++) {
            Node *Ni = &NodeSet[i];
            for (j = 1; j <= Dim; j++) {
                Node *Nj = &NodeSet[j];
                if (Ni != Nj &&
                    Ni->Earliest + Ni->ServiceTime + Ni->C[j] > Nj->Latest)
                    Ni->C[j] = M;
            }
        }
    }
    C = WeightType == EXPLICIT ? C_EXPLICIT : C_FUNCTION;
    D = WeightType == EXPLICIT ? D_EXPLICIT : D_FUNCTION;
    if (ProblemType != CVRP && ProblemType != CVRPTW &&
        ProblemType != CTSP &&
        ProblemType != TSP && ProblemType != ATSP) {
        M = INT_MAX / 2 / Precision;
        for (i = Dim + 1; i <= DimensionSaved; i++) {
            for (j = 1; j <= DimensionSaved; j++) {
                if (j == i)
                    continue;
                if (j == MTSPDepot || j > Dim)
                    NodeSet[i].C[j] = NodeSet[MTSPDepot].C[j] = M;
                NodeSet[i].C[j] = NodeSet[MTSPDepot].C[j];
                NodeSet[j].C[i] = NodeSet[j].C[MTSPDepot];
            }
        }
        if (ProblemType == CCVRP || ProblemType == OVRP)
            for (i = 1; i <= Dim; i++)
                NodeSet[i].C[MTSPDepot] = 0;
    }
    if (Precision > 1 && CostMatrix) {
        for (i = 2; i <= Dim; i++) {
            Node *N = &NodeSet[i];
            for (j = 1; j < i; j++)
                if (N->C[j] * Precision / Precision != N->C[j])
                    eprintf("PRECISION (= %d) is too large", Precision);
        }
    }
    if (SubsequentMoveType == 0) {
        SubsequentMoveType = MoveType;
        SubsequentMoveTypeSpecial = MoveTypeSpecial;
    }
    K = MoveType >= SubsequentMoveType || !SubsequentPatching ?
        MoveType : SubsequentMoveType;
    if (PatchingC > K)
        PatchingC = K;
    if (PatchingA > 1 && PatchingA >= PatchingC)
        PatchingA = PatchingC > 2 ? PatchingC - 1 : 1;
    if (NonsequentialMoveType == -1 ||
        NonsequentialMoveType > K + PatchingC + PatchingA - 1)
        NonsequentialMoveType = K + PatchingC + PatchingA - 1;
    if (PatchingC >= 1 && NonsequentialMoveType >= 4) {
        BestMove = BestSubsequentMove = BestKOptMove;
        if (!SubsequentPatching && SubsequentMoveType <= 5) {
            MoveFunction BestOptMove[] =
            { 0, 0, Best2OptMove, Best3OptMove,
                Best4OptMove, Best5OptMove
            };
            BestSubsequentMove = BestOptMove[SubsequentMoveType];
        }
    } else {
        MoveFunction BestOptMove[] = { 0, 0, Best2OptMove, Best3OptMove,
            Best4OptMove, Best5OptMove
        };
        BestMove = MoveType <= 5 ? BestOptMove[MoveType] : BestKOptMove;
        BestSubsequentMove = SubsequentMoveType <= 5 ?
            BestOptMove[SubsequentMoveType] : BestKOptMove;
    }
    if (MoveTypeSpecial)
        BestMove = BestSpecialOptMove;
    if (SubsequentMoveTypeSpecial)
        BestSubsequentMove = BestSpecialOptMove;
    if (ProblemType == HCP || ProblemType == HPP)
        MaxCandidates = 0;
    if (TraceLevel >= 1) {
        printff("done\n");
        PrintParameters();
    }
    free(LastLine);
    LastLine = 0;
}

void OutputTourResult(LKHOutput* lkhOutput){
    int i, j, n, Forwards;

    n = ProblemType != ATSP ? Dimension : Dimension / 2;
    for (i = 1; i < n && BestTour[i] != 1; i++);
    Forwards = ProblemType == ATSP ||
        BestTour[i < n ? i + 1 : 1] < BestTour[i > 1 ? i - 1 : Dimension];
    for (j = 1; j <= n; j++) {
        lkhOutput->tourResult[j - 1] = BestTour[i];
        if (Forwards) {
            if (++i > n)
                i = 1;
        } else if (--i < 1)
            i = n;
    }
    lkhOutput->tourCost = BestCost;
}

void ReSetLKH(){
    // [TODO 1: 释放 LKH内核 申请(malloc等)的内存, 可以对照源文件 AllocateStructures.c 来做]
    /* tips: 如果在程序结束前没有调用 free，申请的内存将在程序退出时被操作系统回收，
        但这并不算是良好的编程习惯。确保释放内存是管理动态内存的好方法。 */

    // [TODO 2: 重置 LKH内核 的全局变量, 可以对照源文件 LKH.c 来做]
    /* tips: 因为 LKH内核 的运行依赖于全局变量，如果在一个进程内存在多次调用 LKH内核 的话，
        除第一次调用外，后续调用需要重置全局变量，这样才能确保不会对后续调用产生影响 */

    // [TODO 3: 重置 LKH内核 涉及的静态变量, 这些静态变量分布于多个源文件中，需要逐个确认]
    /* tips: 因为 LKH内核 的运行同样依赖于这些静态变量，如果在一个进程内存在多次调用 LKH内核 的话，
        除第一次调用外，后续调用需要重置这些静态变量，这样才能确保不会对后续调用产生影响 */
    /* method: 由于静态变量的作用域仅限于其定义的源文件，无法在其他源文件中直接重新初始化。
        因此，需要使用函数接口：即在定义静态变量的源文件内部，创建一个函数来重新初始化其涉及的静态变量，
        其他源文件就可以通过调用这个函数来实现重新初始化。 */
}

int solveTSP(LKHInput* lkhInput, LKHOutput* lkhOutput){
    /* 参数设定 */
    ReadParameters(lkhInput);
    /* 问题信息设定、数据读入 */
    ReadProblem(lkhInput);
    /* 调用 LKH内核 */
    int ret = LKHmain();
    /* 提取求解结果 */
    OutputTourResult(lkhOutput);
    /* 重置 LKH内核 状态*/
    ReSetLKH();
    return ret;
}

static char *Copy(char *S)
{
    char *Buffer;

    if (!S || strlen(S) == 0)
        return 0;
    Buffer = (char *) malloc(strlen(S) + 1);
    strcpy(Buffer, S);
    return Buffer;
}
static void CheckSpecificationPart()
{
    if (ProblemType == -1)
        eprintf("TYPE is missing");
    if (Dimension < 3)
        eprintf("DIMENSION < 3 or not specified");
    if (WeightType == -1 && !Asymmetric && ProblemType != HCP &&
        ProblemType != HPP && !EdgeWeightType && ProblemType != STTSP)
        eprintf("EDGE_WEIGHT_TYPE is missing");
    if (WeightType == EXPLICIT && WeightFormat == -1 && !EdgeWeightFormat)
        eprintf("EDGE_WEIGHT_FORMAT is missing");
    if (WeightType == EXPLICIT && WeightFormat == FUNCTION)
        eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if (WeightType != EXPLICIT &&
        (WeightType != SPECIAL || CoordType != NO_COORDS) &&
        WeightType != -1 && WeightFormat != -1 && WeightFormat != FUNCTION)
        eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if ((ProblemType == ATSP || ProblemType == SOP) &&
        WeightType != EXPLICIT && WeightType != -1)
        eprintf("Conflicting TYPE and EDGE_WEIGHT_TYPE");
    if (CandidateSetType == DELAUNAY && !TwoDWeightType() &&
        MaxCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = DELAUNAY");
    if (CandidateSetType == NN && !TwoDWeightType() &&
        MaxCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for "
             "CANDIDATE_SET_TYPE = NEAREST-NEIGHBOR");
    if (CandidateSetType == QUADRANT && !TwoDWeightType() &&
        !ThreeDWeightType() && MaxCandidates + ExtraCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = QUADRANT");
    if (ExtraCandidateSetType == QUADRANT && !TwoDWeightType() &&
        !ThreeDWeightType() && ExtraCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for EXTRA_CANDIDATE_SET_TYPE = "
             "QUADRANT");
    if (InitialTourAlgorithm == QUICK_BORUVKA && !TwoDWeightType() &&
        !ThreeDWeightType())
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "QUICK-BORUVKA");
    if (InitialTourAlgorithm == SIERPINSKI && !TwoDWeightType())
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "SIERPINSKI");
    if (DelaunayPartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for DELAUNAY specification");
    if (KarpPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for KARP specification");
    if (KCenterPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for K-CENTER specification");
    if (KMeansPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for K-MEANS specification");
    if (MoorePartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for MOORE specification");
    if (RohePartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for ROHE specification");
    if (SierpinskiPartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for SIERPINSKI specification");
    if (SubproblemBorders && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for BORDERS specification");
    if (InitialTourAlgorithm == MTSP_ALG && Asymmetric)
        eprintf("INTIAL_TOUR_ALGORITHM = MTSP is not applicable for "
                "asymetric problems");
}
static void CreateNodes()
{
    Node *Prev = 0, *N = 0;
    int i;
    if (Dimension <= 0)
        eprintf("DIMENSION is not positive (or not specified)");
    if (Asymmetric) {
        Dim = DimensionSaved;
        DimensionSaved = Dimension + Salesmen - 1;
        Dimension = 2 * DimensionSaved;
    } else if (ProblemType == HPP) {
        Dimension++;
        if (Dimension > MaxMatrixDimension)
            eprintf("DIMENSION too large in HPP problem");
    }
    NodeSet = (Node *) calloc(Dimension + 1, sizeof(Node));
    for (i = 1; i <= Dimension; i++, Prev = N) {
        N = &NodeSet[i];
        if (i == 1)
            FirstNode = N;
        else
            Link(Prev, N);
        N->Id = N->OriginalId = i;
        if (MergeTourFiles >= 1)
            N->MergeSuc = (Node **) calloc(MergeTourFiles, sizeof(Node *));
        N->Earliest = 0;
        N->Latest = INT_MAX;
    }
    Link(N, FirstNode);
}
static int FixEdge(Node * Na, Node * Nb)
{
    if (!Na->FixedTo1 || Na->FixedTo1 == Nb)
        Na->FixedTo1 = Nb;
    else if (!Na->FixedTo2 || Na->FixedTo2 == Nb)
        Na->FixedTo2 = Nb;
    else
        return 0;
    if (!Nb->FixedTo1 || Nb->FixedTo1 == Na)
        Nb->FixedTo1 = Na;
    else if (!Nb->FixedTo2 || Nb->FixedTo1 == Na)
        Nb->FixedTo2 = Na;
    else
        return 0;
    return 1;
}
static void Read_FIXED_EDGES_SECTION(LKHInput* lkhInput){
    Node *Ni, *Nj, *N, *NPrev = 0, *NNext;
    int i, j, Count = 0;

    CheckSpecificationPart();
    if (!FirstNode)
        CreateNodes();
    if (ProblemType == HPP)
        Dimension--;
    int idx = 0;
    while (idx < lkhInput->fixEdgeLen * 2) {
        i = lkhInput->fixEdge[idx++];
        if (i <= 0 || i > (Asymmetric ? Dimension / 2 : Dimension))
            eprintf("FIXED_EDGES_SECTION: Node number out of range: %d",
                    i);
        j = lkhInput->fixEdge[idx++];
        if (j <= 0 || j > (Asymmetric ? Dimension / 2 : Dimension))
            eprintf("FIXED_EDGES_SECTION: Node number out of range: %d",
                    j);
        if (i == j)
            eprintf("FIXED_EDGES_SECTION: Illegal edge: %d to %d", i, j);
        Ni = &NodeSet[i];
        Nj = &NodeSet[Asymmetric ? j + Dimension / 2 : j];
        if (!FixEdge(Ni, Nj))
            eprintf("FIXED_EDGES_SECTION: Illegal fix: %d to %d", i, j);
        /* Cycle check */
        N = Ni;
        Count = 0;
        do {
            NNext = N->FixedTo1 != NPrev ? N->FixedTo1 : N->FixedTo2;
            NPrev = N;
            Count++;
        } while ((N = NNext) && N != Ni);
        if (N == Ni && Count != Dimension)
            eprintf("FIXED_EDGES_SECTION: Illegal fix: %d to %d", i, j);
    }
    if (ProblemType == HPP)
        Dimension++;
}
static void Read_TOUR_SECTION(LKHInput* lkhInput) {
    Node* First = 0, * Last = 0, * N, * Na;
    int i, k;
    
    if (TraceLevel >= 1) {
        printff("Reading TOUR...");
    }
    if (!FirstNode)
        CreateNodes();
    N = FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != FirstNode);
    if (ProblemType == HPP)
        Dimension--;
    if (Asymmetric)
        Dimension = DimensionSaved;
    int idx = 0;
    i = lkhInput->intialTour[idx];
    for (k = 0; k <=Dimension ; k++) {
        if (i <= 0 || i > Dimension)
            eprintf("(TOUR_SECTION) Node number out of range: %d",i);
        N = &NodeSet[i];
        if (N->V == 1 && k != Dimension)
            eprintf("(TOUR_SECTION) Node number occurs twice: %d",N->Id);
        N->V = 1;
        if (k == 0)
            First = Last = N;
        else {
            if (Asymmetric) {
                Na = N + Dimension;
                Na->V = 1;
            }
            else
                Na = 0;
            /* 划分子问题需要 */
            if (!Na)
                (Last->SubproblemSuc = N)->SubproblemPred = Last;
            else {
                (Last->SubproblemSuc = Na)->SubproblemPred = Last;
                (Na->SubproblemSuc = N)->SubproblemPred = Na;
            }
            /* 初始解需要 */
            if (!Na)
                Last->InitialSuc = N;
            else {
                Last->InitialSuc = Na;
                Na->InitialSuc = N;
            }

            Last = N;
        }
        if (k < Dimension) {
            i = lkhInput->intialTour[++idx];
        }
        if (k == Dimension - 1)
            i = First->Id;
    }
  
    N = FirstNode;
    do {    
        if (!N->V)
            printf("TOUR_SECTION: Node is missing: %d", N->Id);
    } while ((N = N->Suc) != FirstNode);
        do {
            if (N->FixedTo1 &&
                N->SubproblemPred != N->FixedTo1
                && N->SubproblemSuc != N->FixedTo1)
                eprintf("Fixed edge (%d, %d) "
                    "does not belong to subproblem tour", N->Id,
                    N->FixedTo1->Id);
            if (N->FixedTo2 && N->SubproblemPred != N->FixedTo2
                && N->SubproblemSuc != N->FixedTo2)
                eprintf("Fixed edge (%d, %d) "
                    "does not belong to subproblem tour", N->Id,
                    N->FixedTo2->Id);
        } while ((N = N->Suc) != FirstNode);
    if (ProblemType == HPP)
        Dimension++;
    if (Asymmetric)
        Dimension *= 2;
    if (TraceLevel >= 1)
        printff("done\n");
}
static void Read_EDGE_WEIGHT_SECTION(LKHInput* lkhInput)
{
    Node *Ni;
    int i, j, n, W;
    double w;
    CheckSpecificationPart();
    if (!FirstNode)
        CreateNodes(); // 创建结点数组，是一个头尾相连的双链表数据结构，FirstNode指向第一个结点
    if (!Asymmetric) {
        CostMatrix = (int *) calloc((size_t) Dimension * (Dimension - 1) / 2,
                                    sizeof(int));
        Ni = FirstNode->Suc;
        do {
            Ni->C =
                &CostMatrix[(size_t) (Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
        }
        while ((Ni = Ni->Suc) != FirstNode);
    } else {
        n = Dimension / 2;
        CostMatrix = (int *) calloc((size_t) n * n, sizeof(int));
        for (Ni = FirstNode; Ni->Id <= n; Ni = Ni->Suc)
            Ni->C = &CostMatrix[(size_t) (Ni->Id - 1) * n] - 1;
    }
    if (ProblemType == HPP)
        Dimension--;
    if (Scale < 1)
        Scale = 1;
    // 写入矩阵中
    for (i = 1; i <= Dim; i++) {
        Ni = &NodeSet[i];
        for (j = 1; j <= Dim; j++) {
            W = lkhInput->adjMat[i - 1][j - 1]; // 对应位置
            if (Asymmetric) {
                Ni->C[j] = W;
                if (j != i && W > M)
                    M = W;
            } else if (j < i)
                Ni->C[j] = W;
        }
    }
    if (ProblemType == HPP)
        Dimension++;
    if (Asymmetric) {
        for (i = 1; i <= DimensionSaved; i++)
            FixEdge(&NodeSet[i], &NodeSet[i + DimensionSaved]);
        Distance = Distance_ATSP;
        WeightType = -1;
    }
}
static int TwoDWeightType()
{
    if (Asymmetric)
        return 0;
    return WeightType == EUC_2D || WeightType == MAX_2D ||
        WeightType == MAN_2D || WeightType == CEIL_2D ||
        WeightType == FLOOR_2D ||
        WeightType == GEO || WeightType == GEOM ||
        WeightType == GEO_MEEUS || WeightType == GEOM_MEEUS ||
        WeightType == ATT || WeightType == TOR_2D ||
        (WeightType == SPECIAL && CoordType == TWOD_COORDS);
}
static int ThreeDWeightType()
{
    if (Asymmetric)
        return 0;
    return WeightType == EUC_3D || WeightType == MAX_3D ||
        WeightType == MAN_3D || WeightType == CEIL_3D ||
        WeightType == FLOOR_3D || WeightType == TOR_3D ||
        WeightType == XRAY1 || WeightType == XRAY2 ||
        (WeightType == SPECIAL && CoordType == THREED_COORDS);
}
