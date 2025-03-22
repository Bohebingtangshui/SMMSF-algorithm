function [IDX]=TwoRoundMerge(clusters,MST_1,MST_3,data,K,clustersCenters)
CenterData=data(clustersCenters,:);
N=size(data,1);
num_C=numel(clusters);

%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
if numel(clusters)>K
clusterSizes = cellfun(@numel, clusters);
thresholdN=N^0.5;
ConnectionShip=GetConnectionShip(clusters,MST_1);%查询所有小簇对在MST_1中是否的邻接簇
CenterDis=pdist2(CenterData,CenterData);         %计算中心点对之间距离
maxDisClu=GetMaxCluDis(CenterData,data,clusters);%计算每个小簇内中心点到最远点距离
MergeFilter=GetMergeFilter(maxDisClu,CenterDis); %计算小簇对是否满足第三个条件
visitedC=zeros(num_C,1);
%遍历簇对，满足三个条件的标记为同一个簇
for i=1:num_C-1
    CurrentLabel=max(visitedC)+1;
    if visitedC(i,1)>0 
        CurrentLabel=visitedC(i,1); 
    else
        visitedC(i,1)=CurrentLabel;
    end
    for j=i+1:num_C
        if ConnectionShip(i,j)>0 && MergeFilter(i,j)>0 && visitedC(j,1)==0 && (clusterSizes(i,1)<thresholdN && clusterSizes(j,1)<thresholdN)
            visitedC(j,1)=CurrentLabel;
        end
    end
end
if visitedC(end,1)==0
    visitedC(end,1)=max(visitedC)+1;
end
%合并
for i=1:max(visitedC)
    clusterIndices = visitedC == i;
    mergedCluster{i,1} = vertcat(clusters{clusterIndices,1});
end
clusters=mergedCluster;
end

%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STEP2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
num_C=numel(clusters);
while(num_C>K)
    InterConnetionRatio=zeros(num_C);
    for i=1:num_C
        for j=(i+1):num_C
            EdgeC1C2=MST_3(clusters{i,1},clusters{j,1});
            [row,~]=find(EdgeC1C2>0);
            if isempty(row)
                continue;
            else
                [VertRatio,d]=computeVertRatio(clusters{i,1},clusters{j,1},MST_3);%计算IC值
            end
            InterConnetionRatio(i,j)=VertRatio*d;
        end
    end
    %寻找最大IC合并,删去空白簇
    maxIC=max(max(InterConnetionRatio));
    [rowIndex,colIndex]=find(InterConnetionRatio==maxIC);
    rowIndex=rowIndex(1);
    colIndex=colIndex(1);
    clusters{rowIndex,1}=[clusters{rowIndex,1};clusters{colIndex,1}];
    clusters{colIndex,1}=[];
    clusters(cellfun(@isempty,clusters))=[];
    num_C=numel(clusters);
end
for i=1:K
    IDX(clusters{i,1})=i;
end
IDX=IDX';
end

%%
%计算IC值
function [VertRatio,d]=computeVertRatio(C1,C2,MST_3)
VertRatio=0;
num1=numel(C1);
num2=numel(C2);
EdgeC1C2=MST_3(C1,C2);%割边集合
[row,col]=find(EdgeC1C2>0);
if isempty(row)
    return;
end
V1=unique(row);%割点1
V2=unique(col);%割点2
numV1=numel(V1);%子簇1数目
numV2=numel(V2);%子簇2数目
VertRatio = (numV1 + numV2 ) / (num1+num2);%第一部分
AvgEdge=sum(EdgeC1C2,'all')/numel(row);%平均割边
V1pt=C1(V1);
V2pt=C2(V2);
restptC1=setdiff(C1,V1pt);
restptC2=setdiff(C2,V2pt);
conntgraphC1=MST_3(V1pt,restptC1);%割点与非割点的边
conntgraphC2=MST_3(V2pt,restptC2);
AvgConnC1=mean(conntgraphC1(conntgraphC1>0));
AvgConnC2=mean(conntgraphC2(conntgraphC2>0));
if AvgConnC1>=AvgEdge
    d1=AvgEdge/AvgConnC1;
else
    d1=AvgConnC1/AvgEdge;
end
if AvgConnC2>=AvgEdge
    d2=AvgEdge/AvgConnC2;
else
    d2=AvgConnC2/AvgEdge;
end
d=min(d1,d2);%第二部分
end

%%
%查询簇对之间是否满足第三个条件
function MergeFilter=GetMergeFilter(maxDisClu,CenterDis)
numClu=numel(maxDisClu);
MergeFilter=zeros(numClu);
for i=1:numClu-1
    for j=i+1:numClu
        if maxDisClu(i,1)+maxDisClu(j,1)>=1.5*CenterDis(i,j)
            MergeFilter(i,j)=1;
        end
    end
end
end
%%
%计算每个小簇内中心点到最远点距离
function [maxDisClu]=GetMaxCluDis(CenterData,data,clusters)
num_Clu=numel(clusters);
maxDisClu = []; % 初始化存储最大距离的数组
for i = 1:num_Clu
    cluster_points = data(clusters{i,1}, :);
    center_point = CenterData(i, :);
    distances_to_center = pdist2(cluster_points, center_point);
    max_distance = max(distances_to_center);
    maxDisClu(end+1, 1) = max_distance;
end
end
%%
%查询所有小簇对在MST_1中是否的邻接簇
function [ConnectionShip]=GetConnectionShip(clusters,MST_1)
ConnectionShip=zeros(size(clusters,1));
for i = 1:size(clusters,1)-1
    for j=i+1:size(clusters,1)
        clu_X=clusters{i,1};
        clu_Y=clusters{j,1};
        MST_XY=MST_1(clu_X,clu_Y);
        positive_elements = MST_XY(MST_XY > 0);
        mean_value = mean(positive_elements);
        if mean_value>0
            ConnectionShip(i,j)=1;
        end
    end
end
end
