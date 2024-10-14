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

static char *Copy(char *S);
static void CreateNodes(void);
static int FixEdge(Node * Na, Node * Nb);
static void Read_FIXED_EDGES_SECTION(LKHInput* lkhInput);
static void Read_TOUR_SECTION(LKHInput* lkhInput);
static void Read_EDGE_WEIGHT_SECTION(LKHInput* lkhInput);

void loadDefaultParam(LKHParameters* p){
    p->TimeSpan = 1; 
    p->ScheduleScoreInSecond = 1000;

    p->ProblemFileName = 0;
    p->PiFileName = 0;
    p->InputTourFileName = 0;
    p->OutputTourFileName = 0;
    p->TourFileName = 0;
    p->CandidateFiles = 0;
    p->MergeTourFiles = 0;

    p->AscentCandidates = 25;//50
    p->BackboneTrials = 0;
    p->Backtracking = 0;
    p->CandidateSetSymmetric = 0;
    p->CandidateSetType = ALPHA;
    p->Crossover = ERXT;
    p->DelaunayPartitioning = 0;
    p->DelaunayPure = 0;
    p->Excess = -1;
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
    p->MoveType = 3;
    p->NonsequentialMoveType = -1;
    p->Optimum = MINUS_INFINITY;
    p->PatchingA = 1;
    p->PatchingC = 0;
    p->PatchingAExtended = 0;
    p->PatchingARestricted = 0;
    p->PatchingCExtended = 0;
    p->PatchingCRestricted = 0;
    p->Precision = 100;
    p->POPMUSIC_InitialTour = 0;
    p->POPMUSIC_MaxNeighbors = 5;
    p->POPMUSIC_SampleSize = 15;
    p->POPMUSIC_Solutions = 15;
    p->POPMUSIC_Trials = 3;
    p->Recombination = IPT;
    p->RestrictedSearch = 1;
    p->RohePartitioning = 0;
    p->Runs = 0;
    p->Seed = 1;
    p->SierpinskiPartitioning = 0;
    p->StopAtOptimum = 1;
    p->Subgradient = 1;
    p->SubproblemBorders = 0;
    p->SubproblemsCompressed = 0;
    p->SubproblemSize = 0;
    p->SubsequentMoveType = 0;
    p->SubsequentPatching = 1;
    p->TimeLimit = DBL_MAX;
    p->TotalTimeLimit = DBL_MAX;
    p->TraceLevel = 1;
    p->MaxMatrixDimension = 20006;
}

void ReadParameters(LKHInput* lkhInput)
{
     LKHParameters* lkhParam = lkhInput->lkhParameters;
    // 从结构体获取参数
    TimeSpan = lkhParam->TimeSpan; 
    ScheduleScoreInSecond = lkhParam->ScheduleScoreInSecond;

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
    CandidateSetSymmetric = lkhParam->CandidateSetSymmetric;
    CandidateSetType = lkhParam->CandidateSetType;
    Crossover = lkhParam->Crossover;
    DelaunayPartitioning = lkhParam->DelaunayPartitioning;
    DelaunayPure = lkhParam->DelaunayPure;
    Excess = lkhParam->Excess;
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
    NonsequentialMoveType = lkhParam->NonsequentialMoveType;
    Optimum  = lkhParam->Optimum ;
    PatchingA = lkhParam->PatchingA;
    PatchingC = lkhParam->PatchingC;
    PatchingAExtended = lkhParam->PatchingAExtended;
    PatchingARestricted = lkhParam->PatchingARestricted;
    PatchingCExtended = lkhParam->PatchingCExtended;
    PatchingCRestricted = lkhParam->PatchingCRestricted;
    Precision = lkhParam->Precision;
    POPMUSIC_InitialTour = lkhParam->POPMUSIC_InitialTour;
    POPMUSIC_MaxNeighbors = lkhParam->POPMUSIC_MaxNeighbors;
    POPMUSIC_SampleSize = lkhParam->POPMUSIC_SampleSize;
    POPMUSIC_Solutions = lkhParam->POPMUSIC_Solutions;
    POPMUSIC_Trials = lkhParam->POPMUSIC_Trials;
    Recombination = lkhParam->Recombination;
    RestrictedSearch = lkhParam->RestrictedSearch;
    RohePartitioning = lkhParam->RohePartitioning;
    Runs = lkhParam->Runs;
    Seed = lkhParam->Seed;
    SierpinskiPartitioning = lkhParam->SierpinskiPartitioning;
    StopAtOptimum = lkhParam->StopAtOptimum;
    Subgradient = lkhParam->Subgradient;
    SubproblemBorders = lkhParam->SubproblemBorders;
    SubproblemsCompressed = lkhParam->SubproblemsCompressed;
    SubproblemSize = lkhParam->SubproblemSize;
    SubsequentMoveType = lkhParam->SubsequentMoveType;
    SubsequentPatching = lkhParam->SubsequentPatching;
    TimeLimit = lkhParam->TimeLimit;
    TotalTimeLimit = lkhParam->TotalTimeLimit;
    TraceLevel = lkhParam->TraceLevel;
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
    /* NAME */
    Name = Copy("IOSchedule");
    /* TYPE */
    Type = Copy("ATSP");
    ProblemType = ATSP;
    /* DIMENSION */
    Dimension = lkhInput->matDimension;;
    DimensionSaved = Dimension;
    /* EDGE_WEIGHT_TYPE */
    EdgeWeightType = Copy("EXPLICIT");
    WeightType = EXPLICIT;
    Distance = Distance_EXPLICIT;
    /* EDGE_WEIGHT_FORMAT */
    EdgeWeightFormat = Copy("FULL_MATRIX");
    WeightFormat = FULL_MATRIX;
    /* EDGE_DATA_FORMAT */
    EdgeDataFormat = NULL; 
    /* NODE_COORD_TYPE */
    NodeCoordType = Copy("NO_COORDS");
    CoordType = NO_COORDS;
    /* DISPLAY_DATA_TYPE */
    DisplayDataType = Copy("NO_DISPLAY");
    /* GRID_SIZE */
    GridSize = 1000000.0; // no meaning in here 

    /* (2) The data part 根据上方定义，进行实际数据输入 */
    /* NODE_COORD_SECTION */
    // NULL
    /* EDGE_DATA_SECTION */
    // NULL
    /* FIXED_EDGES_SECTION */
    if(lkhInput->fixEdgeLen != 0)
        Read_FIXED_EDGES_SECTION(lkhInput);
    /* DISPLAY_DATA_SECTION */
    // NULL
    /* TOUR_SECTION */
    if(lkhInput->intialTour != 0)
        Read_TOUR_SECTION(lkhInput);
    /* EDGE_WEIGHT_SECTION */
    Read_EDGE_WEIGHT_SECTION(lkhInput);

    Swaps = 0;

    /* Adjust parameters 调整参数。比如某个参数的默认值依赖于问题类型或输入数据，
       在此前(ReadParameter 函数中)无法确定值是多少，仅是设定为一个默认值标记，
       所以在此处(ReadProblem 函数最后)其默认值才被实际的设定 */
    int i, K;
    if (Seed == 0)
        Seed = (unsigned) (time(0) * (size_t) (&Seed));
    if (Precision == 0)
        Precision = 100;
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
            Excess = 1.0 / Dimension;
        if (MaxTrials == -1)
            MaxTrials = Dimension;
        HeapMake(Dimension);
    }
    if (POPMUSIC_MaxNeighbors > Dimension - 1)
        POPMUSIC_MaxNeighbors = Dimension - 1;
    if (POPMUSIC_SampleSize > Dimension)
        POPMUSIC_SampleSize = Dimension;
   
    // if (Precision > 1 && (WeightType == EXPLICIT || ProblemType == ATSP)) {
    //     int j, n = ProblemType == ATSP ? Dimension / 2 : Dimension;
    //     for (i = 2; i <= n; i++) {
    //         Node *N = &NodeSet[i];
    //         for (j = 1; j < i; j++)
    //             if (N->C[j] * Precision / Precision != N->C[j])
    //                 eprintf("PRECISION (= %d) is too large", Precision);
    //     }
    // }
    C = WeightType == EXPLICIT ? C_EXPLICIT : C_FUNCTION;
    D = WeightType == EXPLICIT ? D_EXPLICIT : D_FUNCTION;
    if (SubsequentMoveType == 0)
        SubsequentMoveType = MoveType;
    K = MoveType >= SubsequentMoveType
        || !SubsequentPatching ? MoveType : SubsequentMoveType;
    if (PatchingC > K)
        PatchingC = K;
    if (PatchingA > 1 && PatchingA >= PatchingC)
        PatchingA = PatchingC > 2 ? PatchingC - 1 : 1;
    if (NonsequentialMoveType == -1 ||
        NonsequentialMoveType > K + PatchingC + PatchingA - 1)
        NonsequentialMoveType = K + PatchingC + PatchingA - 1;
    if (PatchingC >= 1) {
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
    if (TraceLevel >= 1) {
        printff("done\n");
        PrintParameters();
    } 
}

void OutputTourResult(LKHOutput* lkhOutput){
    int j, n = DimensionSaved;
    for (j = 1; j <= n; j++) {
        lkhOutput->tourResult[j - 1] = BestTour[j];
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
    int ret = LKHmain(lkhInput->scheduleStartTime);
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
static void CreateNodes()
{
    Node *Prev = 0, *N = 0;
    int i;

    if (Dimension <= 0)
        eprintf("DIMENSION is not positive (or not specified)");
    if (ProblemType == ATSP)
        Dimension *= 2;
    NodeSet = (Node *) calloc(Dimension + 1, sizeof(Node));
    for (i = 1; i <= Dimension; i++, Prev = N) {
        N = &NodeSet[i];
        if (i == 1)
            FirstNode = N;
        else
            Link(Prev, N);
        N->Id = i;
        if (MergeTourFiles >= 1)
            N->MergeSuc = (Node **) calloc(MergeTourFiles, sizeof(Node *));
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

    if (!FirstNode)
        CreateNodes();
    int idx = 0;
    while (idx < lkhInput->fixEdgeLen * 2) {
        i = lkhInput->fixEdge[idx++];
        if (i <= 0 || i > DimensionSaved)
            eprintf("FIXED_EDGES_SECTION: Node number out of range: %d", i);
        j = lkhInput->fixEdge[idx++];
        if (j <= 0 || j > DimensionSaved)
            eprintf("FIXED_EDGES_SECTION: Node number out of range: %d", j);
        if (i == j)
            eprintf("FIXED_EDGES_SECTION: Illegal edge: %d to %d", i, j);
        Ni = &NodeSet[i];
        Nj = &NodeSet[j + DimensionSaved];
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

    int idx = 0;
    i = lkhInput->intialTour[idx];
    for (k = 0; k <= DimensionSaved ; k++) {
        if (i <= 0 || i > DimensionSaved)
            eprintf("(TOUR_SECTION) Node number out of range: %d",i);
        N = &NodeSet[i];
        if (N->V == 1 && k != DimensionSaved)
            eprintf("(TOUR_SECTION) Node number occurs twice: %d",N->Id);
        N->V = 1;
        if (k == 0)
            First = Last = N;
        else {
            Na = N + DimensionSaved;
            Na->V = 1;
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
        if (k < DimensionSaved) {
            i = lkhInput->intialTour[++idx];
        }
        if (k == DimensionSaved - 1)
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
    if (TraceLevel >= 1)
        printff("done\n");
}
static void Read_EDGE_WEIGHT_SECTION(LKHInput* lkhInput)
{
    Node *Ni, *Nj;
    int i, j, n, W;

    if (!FirstNode)
        CreateNodes(); // 创建结点数组，是一个头尾相连的双链表数据结构，FirstNode指向第一个结点

    /* ProblemType == ATSP && WeightFormat == FULL_MATRIX */
    // 为矩阵申请空间
    n = DimensionSaved;
    CostMatrix = (int *) calloc((size_t) n * n, sizeof(int));
    for (Ni = FirstNode; Ni->Id <= n; Ni = Ni->Suc)
        Ni->C = &CostMatrix[(size_t) (Ni->Id - 1) * n] - 1;
    // 写入矩阵中
    for (i = 1; i <= n; i++) {
        Ni = &NodeSet[i];
        for (j = 1; j <= n; j++) {
            W = lkhInput->adjMat[i - 1][j - 1]; // 对应位置
            if (j != i && W > INT_MAX / 2 / Precision)
                eprintf("EDGE_WEIGHT_SECTION: "
                        "Weight %d > INT_MAX / 2 / PRECISION", W);
            Ni->C[j] = W;
            if (i != j && W > M)
                M = W;
        }
        Nj = &NodeSet[i + n];
        FixEdge(Ni, Nj);
    }
    Distance = Distance_ATSP;
    WeightType = -1;
}
