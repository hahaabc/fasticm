include("graph.jl")
function init(G,T,p,len,del)
    n=G.n;
    m=G.m;
    t=size(T)[1];
    wlks=0;
    wlkedge=0;
    wend=0;
    wuse=0;
    wgood=0;
    wleng=0;
    wall=0;
    W=[];
    We=[];
    xyinT=zeros(Int,n)
    xyinT[T].=1;
    for j=1:m
        if del[j]==1
            continue
        end
        # p=rand([1,2])
        x=G.u[j];
        y=G.v[j];
        for i=1:p
            w1=Int64[];
            w2=Int64[];
            xx=x;yy=y;
            lx=1;ly=1;
            push!(w1,xx);push!(w2,yy);
            for ll in 1:len
                if xyinT[xx]==1
                    lx=ll;
                    break
                end
                if length(G.nbr[xx])==0
                    println(xx)
                end
                nx = rand(G.nbr[xx])
                push!(w1,nx);
                xx=nx;
            end
            for ll in 1:len
                if xyinT[yy]==1
                    ly=ll;
                    break
                end
                ny = rand(G.nbr[yy])
                push!(w2,ny);
                yy=ny;
            end
            ii=w1[lx];jj=w2[ly];

            wlks+=1;
            if xyinT[ii]*xyinT[jj]==1
                # adj[ii,jj]+=1/(lx+ly-1)
                # adj[jj,ii]+=1/(lx+ly-1)
                wuse+=1;
                reverse!(w2)
                append!(w2,w1)
                push!(W,w2);
                push!(We,ly);
                wall+=lx+ly-1;
            end
        end
    end

    return W,We,wuse/wlks,wall/wuse;
end
