using StaticArrays

mi=0;ma=0;
t=zeros(4,4);
J1=ones(4,4)/4;
J2=ones(3,3)/3;
function choose_node(G,L)
    n=G.n
    J1=ones(n,n)/n;
    Ld = inv(L+J1)-J1;#mppinv(L)
    z=zeros(n);
    for i=1:G.n
        z[i]=Ld[i,i];
    end
    ans=argmax(z)
    return ans
end
function adj2lap(adj::Array{Float64,2})
    s=size(adj)[1]
    for i=1:s
        adj[i,i]=0
    end
    ds=sum(adj,dims=1)
    for i=1:s
        adj[i,i]=ds[i]
        for j=i+1:s
            adj[i,j]=-adj[i,j]
            adj[j,i]=-adj[j,i]
        end
    end
end
function getMTXER(l::Array{Float64,2})
    adj2lap(l)
    s=size(l)[1];
    # println(l)
    invl=inv(l+J1)-J1
    er=invl[1,1]+invl[s,s]-2*invl[1,s]
    return er
end
function get3MTXER(l::Array{Float64,2})
    adj2lap(l)
    s=size(l)[1];
    invl=inv(l+J2)-J2
    er=invl[1,1]+invl[s,s]-2*invl[1,s]
    return er
end
function getlittleER(ll::Array{Float64,2},oll::Array{Float64,2})
    # println()
    # println(ll)
    # println(oll)
    s=size(ll)[1];
    for i=1:s
        for j=i+1:s
            if ll[i,j]<1e-8
                ll[i,j]=0
            end
        end
    end
    nzr=0;
    if s==3
        ER=0;
        if ll[2,3]!=0&&ll[1,2]!=0&&ll[1,3]!=0
            ER=1/(1/(1/ll[1,2]+1/ll[2,3])+ll[1,3])
        elseif ll[2,3]!=0&&ll[1,2]!=0&&ll[1,3]==0
            ER=1/ll[1,2]+1/ll[2,3]
        elseif ll[2,3]!=0&&ll[1,2]==0&&ll[1,3]!=0
            ER=1/(ll[1,3])
        elseif ll[2,3]==0&&ll[1,2]!=0&&ll[1,3]!=0
            ER=1/(ll[1,3])
        end
        # ER=get3MTXER(oll)
    elseif s==4
        # t[:,:]=ll[:,:]
        # t[t.>0]=1 ./t[t.>0]
        for i=1:4
            for j=i+1:4
                # t[i,j]=ll[i,j]
                if ll[i,j]>1e-8
                    nzr+=1;
                    ll[i,j]=1/ll[i,j]
                end
            end
        end
        ER=0;
        # ER=getMTXER(oll)
        if nzr==6||nzr==5
            ER=getMTXER(oll)
            # ER=0;
        elseif nzr==4
            if ll[1,2]==0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=ll[1,4]
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=getMTXER(oll)
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=1/(1/(ll[1,4])+1/(ll[1,3]+ll[3,4]))
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]!=0
                ER=1/(1/(ll[1,4])+1/(ll[1,3]+ll[3,4]))
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]==0
                ER=1/(1/(ll[1,4])+1/(ll[1,3]+ll[2,3]+ll[2,4]))
            elseif ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]==0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=getMTXER(oll)
            elseif ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=1/(1/(ll[1,2]+ll[2,4])+1/(ll[1,4]))
            elseif ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]!=0
                ER=1/(1/(ll[1,2]+ll[2,3]+ll[3,4])+1/(ll[1,4]))
            elseif ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]==0
                ER=1/(1/(ll[1,2]+ll[2,4])+1/(ll[1,4]))
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=1/(1/(ll[1,2]+ll[2,4])+1/(ll[1,3]+ll[3,4]))
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]!=0
                ER=getMTXER(oll)
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]==0
                ER=getMTXER(oll)
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]==0&&ll[2,4]==0&&ll[3,4]!=0
                ER=1/(1/(ll[1,3]+ll[3,4])+1/(ll[1,4]))
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]==0
                ER=1/(1/(ll[1,2]+ll[2,4])+1/(ll[1,4]))
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]==0
                ER=ll[1,4]
            end
        elseif nzr==3
            if ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]==0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=1/(1/(ll[1,2]+ll[2,4]))
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]==0
                ER=1/(1/(ll[1,2]+ll[2,4]))
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]==0&&ll[2,4]==0&&ll[3,4]!=0
                ER=1/(1/(ll[1,3]+ll[3,4]))
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=1/(1/(ll[1,3]+ll[3,4]))
            elseif ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]==0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]==0
                ER=1/(1/(ll[1,2]+ll[2,4]))
            elseif ll[1,2]!=0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]==0&&ll[2,4]==0&&ll[3,4]==0
                ER=ll[1,4]
            elseif ll[1,2]==0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]!=0
                ER=ll[1,4]
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]!=0
                ER=1/(1/(ll[1,3]+ll[3,4]))
            elseif ll[1,2]==0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]==0
                ER=ll[1,4]
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]==0
                ER=ll[1,4]
            elseif ll[1,2]==0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]!=0
                ER=ll[1,4]
            elseif ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]==0
                ER=ll[1,4]
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]==0&&ll[2,3]!=0&&ll[2,4]!=0&&ll[3,4]==0
                ER=1/(1/(ll[1,3]+ll[2,3]+ll[2,4]))
            elseif ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]==0&&ll[2,3]!=0&&ll[2,4]==0&&ll[3,4]!=0
                ER=1/(1/(ll[1,2]+ll[2,3]+ll[3,4]))
            elseif ll[1,2]==0&&ll[1,3]!=0&&ll[1,4]!=0&&ll[2,3]==0&&ll[2,4]!=0&&ll[3,4]==0
                ER=ll[1,4]
            elseif ll[1,2]!=0&&ll[1,3]==0&&ll[1,4]!=0&&ll[2,3]==0&&ll[2,4]==0&&ll[3,4]!=0
                ER=ll[1,4]
            end
        end
    #     # t2=time()
    end
    if ER==0
        # println(s)
        # println(nzr)
        # println(ll)
    end
    return ER
end
function DFS(ll, i,cnt,vis)
    vis[i] = 1;
    cnt[1]+=1;
    for j=1:size(ll)[1]
        if i>j
            mi=j
            ma=i
        else
            mi=i;
            ma=j;
        end
        if ll[mi,ma]>1e-6 && !vis[j]
            DFS(ll, j,cnt,vis);
        end
    end
end
function isconnected(ll)
    s=size(ll)[1];
    cnt=[0];
    vis=falses(s);
    DFS(ll,1,cnt,vis);
    visnum=0;
    return cnt[1]==s
end
function get_init_Conduct(W,We,walk_num,S)
    result=zeros(Float64,n)#conduct to other nodes
    for i=1:walk_num
        global t
        # println(fout,W[i])
        # println(fout,We[i],' ',W[i][We[i]],' ',W[i][We[i]+1])
        # println(W[i])
        cur_walk_len=length(W[i])

        for j=1:We[i]
            vis_l[W[i][j]]=j;
        end
        for j=cur_walk_len:-1:(We[i]+1)
            vis_r[W[i][j]]=j;#绝对位置
        end

        if W[i][cur_walk_len] in S
            for j = 2:We[i]
                t=W[i][j];
                if j==vis_l[t]&&vis_r[t]==0#&&t!=S[1]
                    result[t]+=1/(cur_walk_len-vis_l[t])
                end
            end
            if !(W[i][1] in S)&&1==vis_l[W[i][1]]&&vis_r[W[i][1]]==0
                result[W[i][1]]+=1/(cur_walk_len-1)
            end
        end
        if W[i][1] in S
            for j = We[i]+1:cur_walk_len-1
                t=W[i][j];
                if j==vis_r[t]&&vis_l[t]==0#&&t!=S[1]
                    result[t]+=1/(vis_r[t]-1)
                end
            end
            if !(W[i][cur_walk_len] in S)&&cur_walk_len==vis_r[W[i][cur_walk_len]]&&vis_l[W[i][cur_walk_len]]==0
                result[W[i][1]]+=1/(cur_walk_len-1)
            end
        end
        vis_l[W[i][1:cur_walk_len]].=0;
        vis_r[W[i][1:cur_walk_len]].=0;

    end
    return result
end
function exact_suR()
    #exact computing
    tmp=0;
    J=ones(n,n)/n;
    invL=inv(L+J)-J;
    for j=1:n
        if j!=S[1]
            global tmp+=invL[S[1],S[1]]+invL[j,j]-2*invL[S[1],j];
        end
    end
    return tmp
end

function lap_solver(G,F)
    #laplacian solver
    println("laplacian solver")
    L=lapsp(G);
    f=approxchol_sddm(L[F,F],tol=0);
    M=100;
    tmp=0;
    for i=1:M
        x=randn(length(F))
        tmp+=1/M*(x'*f(x))
    end
    return tmp
end
#二分查找最后一个小于等于某个数,若相等则返回最右边的
function binarySearchFirstLessEqual(vec::Vector{Int64}, target::Int64)
    l = 1;r = length(vec);
	while l <= r
		mid=Int(round(l+(r-l)/2));
        if vec[mid] > target
			r = mid - 1;
		else
			l = mid + 1;
        end
	end
	if r < 0 #如果要找的是最后一个等于target的位置，则条件是 if(r<0 || vec[r]!=target)
		return -1;
    end
	return r;
end
#二分查找第一个大于等于某个数,若相等则返回最左边的
function binarySearchFirstGreaterEqual(vec::Vector{Int64}, target::Int64)
	l = 1;r = length(vec);
	while l <= r
		mid=Int(ceil(l+(r-l)/2));
		if vec[mid] < target
			l = mid + 1;
		else
			r = mid - 1;
		end
    end
	if l > length(vec)    #如果要找的是第一个等于target的位置，则条件是 if(l >= vec.size() || vec[l]!=target)
		return -1;
    end
	return l;
end

function get_real(G,S,del)
    n=G.n;m=G.m
    res=zeros(G.m);
    inir=0
    inirl=[]
    J1=ones(n,n)/n;
    L=lapsp(G);
    #new L
    for i=1:G.m
        if del[i]!=0
        x=G.u[i];y=G.v[i]
        L[x,x]-=1;L[y,y]-=1;L[x,y]=0;L[y,x]=0;
        end
    end
    invL=inv(L+J1)-J1;
    for j=1:n
        if j!=S[1]
            inir+=invL[S[1],S[1]]+invL[j,j]-2*invL[S[1],j];
            push!(inirl,invL[S[1],S[1]]+invL[j,j]-2*invL[S[1],j];)
        else
            push!(inirl,0)
        end
    end

    for i=1:G.m
        if del[i]!=0
            res[i]=0
            continue
        end
        L=lapsp(G);
        for j=1:G.m
            if del[j]!=0
            xx=G.u[j];yy=G.v[j]
            L[xx,xx]-=1;L[yy,yy]-=1;L[xx,yy]=0;L[yy,xx]=0;
            end
        end
        x=G.u[i];y=G.v[i];
        L[x,x]-=1;L[y,y]-=1;L[x,y]=0;L[y,x]=0;
        tmp1=0;
        J1=ones(n,n)/n;
        if !isconnected(-L)
            res[i]=-1e10
            continue
        end
        invL=inv(L+J1)-J1;
        for j=1:n
            if j!=S[1]
                tmp1+=invL[S[1],S[1]]+invL[j,j]-2*invL[S[1],j];
            end
        end
        res[i]=tmp1;
        L[x,x]+=1;L[y,y]+=1;L[x,y]=-1;L[y,x]=-1;
    end
    return inir,res
end
