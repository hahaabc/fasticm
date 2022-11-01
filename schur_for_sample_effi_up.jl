include("graph.jl")
include("core.jl")
include("schur.jl")
include("functions.jl")
include("logw.jl")

using Dates
using LinearAlgebra
fname = open("filename.txt", "r")
str = readline(fname);
nn = parse(Int, str);

for nnn=1:nn
str = readline(fname);
str = split(str);
sclog = open("sclog.txt", "a")
logw(sclog)
logw(sclog, String(fill('*',60)))
logw(sclog, now())
for ro =1:1
for it=1:1
G = get_graph(str[1]);
on = G.n;
om = G.m;
Gc = findconnect(G)
G = Gc;
n = G.n;
m = G.m;

beta = 0.05;
eps = 1;
p = Int(round(log(n) / eps^2));
p = 10;
len = Int(round(1 / beta^2));
len = 200;
selc_num = 1;
kmax = 10;
coeffi = 10;

mode="sc";#"sc,re,bo"
cc = zeros(Int, n)
logw(sclog, str[1], ' ', on, ' ', om, ' ', n, ' ', m)
logw(sclog, "p=", p, ",len=", len, ",selc_num=", selc_num, ",kmax=", kmax);
################Init
L = lapsp(G);
d = zeros(n);
d = [L[i, i] for i = 1:n];
S = Int[];
for i=1:selc_num
    union!(S,rand(1:G.n))
end
G1=G;G2=G;
G=G1
n = G.n;
m = G.m;
del = falses(m)
del2 = falses(m)
L = lapsp(G);
d = zeros(n);
V = union(1:n)
sample_num = coeffi * Int(round(sqrt(n)));
sample_V = rand(V, sample_num)
unique!(sample_V)
sample_num = length(sample_V)
sample_ratio=sample_num/n
logw(sclog, "unionGinfor:", G.n, ',', G.m)
logw(sclog, "sample ", sample_num, " of ", G.n, " nodes, ","sample ratio ",sample_ratio)
logw(sclog, "S:", S[1])
close(sclog)
xyinSam = falses(n)
xyinSam[sample_V] .= true;
F = setdiff(V, S)
deleted=Int64[]
f=open("edges.txt","w")
for i=1:G.m
    u1=G.u[i];u2=G.v[i]
end
close(f)
for e in deleted
    xx = G.u[e];
    yy = G.v[e];
    setdiff!(G.nbr[xx], yy);
    setdiff!(G.nbr[yy], xx);
end
mytt=0;
for k = 1:kmax
    W, We, rat, lavg = init(G, S, p, len, del)
    walk_num = length(We)
    average_len = 0.0
    max_len = 0.0
    min_len = 100000000.0
    for i = 1:walk_num
        average_len += length(W[i])
        if length(W[i]) > max_len
            max_len = length(W[i])
        end
        if length(W[i]) < min_len
            min_len = length(W[i])
        end
    end
    sclog = open("sclog.txt", "a")
    logw(sclog, "")
    logw(sclog, "average:", average_len / walk_num," max:", max_len," min:", min_len," num:", walk_num)
    close(sclog)
    cdct = zeros(n)
    iniRL = zeros(n)
    vis_l = zeros(Int, n)
    vis_r = zeros(Int, n)
    app_ll = trues(n);
    app_rr = trues(n);
    loc_l1 = [Int64[] for i = 1:n]
    loc_l2 = [Int64[] for i = 1:n]
    loc_r1 = [Int64[] for i = 1:n]
    loc_r2 = [Int64[] for i = 1:n]
    lvalue = [Float64[] for i = 1:length(W)]
    rvalue = [Float64[] for i = 1:length(W)]
    lids = [Int64[] for i = 1:length(W)]
    rids = [Int64[] for i = 1:length(W)]
    lns = [Int64[] for i = 1:length(W)]
    rns = [Int64[] for i = 1:length(W)]
    lids1 = [Int64[] for i = 1:length(W)]
    rids1 = [Int64[] for i = 1:length(W)]
    inSamids = [Int64[] for i = 1:length(W)]
    inSamidsLenn = Int64[]

    for i = 1:length(W)
        lenn = size(W[i])[1]
        vis_l[W[i][2:We[i]]] = 2:We[i]
        vis_r[W[i][(lenn-1):-1:(We[i]+1)]] = (lenn-1):-1:(We[i]+1)
        cdct[W[i][2:We[i]]] +=
            (
                (vis_l[W[i][2:We[i]]] .== 2:We[i]) .*
                (vis_r[W[i][2:We[i]]] .== 0) .* (xyinSam[W[i][2:We[i]]])
            ) ./ (lenn-2:-1:lenn-We[i])
        cdct[W[i][lenn-1:-1:We[i]+1]] +=
            (
                (vis_r[W[i][lenn-1:-1:We[i]+1]] .== (lenn-1:-1:We[i]+1)) .*
                (vis_l[W[i][lenn-1:-1:We[i]+1]] .== 0) .*
                xyinSam[W[i][lenn-1:-1:We[i]+1]]
            ) ./ (lenn-2:-1:We[i])
        for j = 2:We[i]
            if vis_l[W[i][j]] == j && xyinSam[W[i][j]]
                push!(loc_l1[W[i][j]], i)
                push!(loc_l2[W[i][j]], j)
                push!(lids1[i], j)
                if vis_r[W[i][j]] == 0
                    push!(lvalue[i], 1 / (lenn - j))
                    push!(lids[i], j)
                    push!(lns[i], W[i][j])
                end
            end
        end
        for j = We[i]+1:lenn-1
            if vis_r[W[i][j]] == j && xyinSam[W[i][j]]
                push!(loc_r1[W[i][j]], i)
                push!(loc_r2[W[i][j]], j)
                push!(rids1[i], j)
                if vis_l[W[i][j]] == 0
                    push!(rvalue[i], 1 / (j - 1))
                    push!(rids[i], j)
                    push!(rns[i], W[i][j])
                end
            end
        end
        vis_l[W[i][1:We[i]]] .= 0
        vis_r[W[i][(We[i]+1):lenn]] .= 0
    end

    for i = 1:n
        push!(loc_l1[i], length(W) + 100)
        push!(loc_l2[i], length(W) + 100)
        push!(loc_r1[i], length(W) + 100)
        push!(loc_r2[i], length(W) + 100)
    end

    cdct[sample_V] ./= p
    for i = 1:n
        if cdct[i] > 0 && xyinSam[i]
            iniRL[i] = 1 / cdct[i]
        end
    end
    sclog = open("sclog.txt", "a")
    logw(sclog,"")
    app_inir = sum(iniRL)/sample_ratio
    logw(sclog, "app_inir: ", app_inir)
    close(sclog)
    res_change = zeros(m)
    weight_upd = zeros(4, 4, n)
    n24 = zeros(Int, n)
    change_id = Int64[]
    tot = zeros(4, 4)
    ll = zeros(4, 4)
    lll = zeros(3, 3)
    oll = zeros(4, 4)
    olll = zeros(3, 3)
    J4 = ones(4, 4) / 4
    J3 = ones(3, 3) / 3
    change_id = falses(n)
    for edge = 1:G.m
        tot .= 0;totij = 0;
        if del[edge] != 0
            res_change[edge]=-10000000
            continue
        end
        x = G.u[edge];y = G.v[edge];
        if y == S[1]
            y = x;x = S[1];
        end
        xsam = false;ysam = false;
        if xyinSam[x]
            xsam = true
        end
        if xyinSam[y]
            ysam = true
        end
        weight_upd[:, :, sample_V] .= 0
        xyinSam[x] = 1;xyinSam[y] = 1;
        il = 1;ir = 1;jl = 1;jr = 1;cur_wid = 0;n24[x] = 2;n24[y] = 3;n24[S] .= 1;
        while min(loc_l1[x][il], loc_r1[x][ir], loc_l1[y][jl], loc_r1[y][jr]) <=
              length(W)
            cur_wid =
                min(loc_l1[x][il], loc_r1[x][ir], loc_l1[y][jl], loc_r1[y][jr])
            lenn = length(W[cur_wid])
            lft = 1
            rht = lenn
            if loc_l1[x][il] == cur_wid && lft < loc_l2[x][il]
                lft=loc_l2[x][il]
            end
            if loc_l1[y][jl] == cur_wid && lft< loc_l2[y][jl]
                lft=loc_l2[y][jl]
            end
            if loc_r1[x][ir] == cur_wid && rht > loc_r2[x][ir]
                rht = loc_r2[x][ir]
            end
            if loc_r1[y][jr] == cur_wid && rht > loc_r2[y][jr]
                rht = loc_r2[y][jr]
            end
            app_ll[W[cur_wid][lft:We[cur_wid]]] .= 0
            app_rr[W[cur_wid][We[cur_wid]+1:rht]] .= 0
            n24l = n24[W[cur_wid][lft]]
            n24r = n24[W[cur_wid][rht]]
            if x != S[1]
                if W[cur_wid][rht] == S[1] &&(W[cur_wid][lft] == x || W[cur_wid][lft] == y)
                    idx = binarySearchFirstLessEqual(lids[cur_wid], lft - 1)
                    weight_upd[1, 4, lns[cur_wid][1:idx]] -=
                        lvalue[cur_wid][1:idx]
                    tt = unique(W[cur_wid][lft+1:lenn-1])
                    v = 1 / (lenn - lft)

                    weight_upd[1, n24l, tt[xyinSam[tt]]] .-= v
                    weight_upd[1, 4, rns[cur_wid]] -= rvalue[cur_wid]

                    bs = rids1[cur_wid][app_ll[W[cur_wid][rids1[cur_wid]]]]
                    weight_upd[n24l, 4, W[cur_wid][bs]] += 1 ./ (bs .- lft)
                    change_id[W[cur_wid][bs]] .= 1
                    if loc_l1[x][il] == loc_l1[y][jl]
                        tmp = loc_l2[x][il]
                        local_node = x
                        if W[cur_wid][lft] == x
                            tmp = loc_l2[y][jl]
                            local_node = y
                        end
                        v = 1 / (lenn - tmp)
                        weight_upd[1, 5-n24l, local_node] += v
                        tot[1, 5-n24l] += v
                    end
                elseif W[cur_wid][lft] == S[1] &&
                       (W[cur_wid][rht] == x || W[cur_wid][rht] == y)
                    idx = binarySearchFirstGreaterEqual(rids[cur_wid], rht + 1)
                    if idx != -1
                        weight_upd[1,4,
                            rns[cur_wid][idx:length(rids[cur_wid])],
                        ] -= rvalue[cur_wid][idx:length(rids[cur_wid])]
                    end
                    tt = unique(W[cur_wid][2:rht-1])
                    v = 1 / (rht - 1)
                    weight_upd[1, n24r, tt[xyinSam[tt]]] .-= v

                    weight_upd[1, 4, lns[cur_wid]] -= lvalue[cur_wid]

                    bs = lids1[cur_wid][app_rr[W[cur_wid][lids1[cur_wid]]]]
                    weight_upd[n24r, 4, W[cur_wid][bs]] += 1 ./ (rht .- bs)
                    change_id[W[cur_wid][bs]] .= 1
                    if loc_r1[x][ir] == loc_r1[y][jr]
                        tmp = loc_r2[x][ir]
                        local_node = x
                        if W[cur_wid][rht] == x
                            tmp = loc_r2[y][jr]
                            local_node = y
                        end
                        v = 1 / (tmp - 1)
                        weight_upd[1, 5-n24r, local_node] += v
                        tot[1, 5-n24r] += v
                    end
                elseif W[cur_wid][lft] == x && W[cur_wid][rht] == x ||
                       W[cur_wid][lft] == y && W[cur_wid][rht] == y#2 x or y
                    idx = binarySearchFirstLessEqual(lids[cur_wid], lft - 1)
                    weight_upd[1, 4, lns[cur_wid][1:idx]] -=
                        lvalue[cur_wid][1:idx]

                    idx = binarySearchFirstGreaterEqual(lids[cur_wid], lft + 1)
                    if idx != -1
                        weight_upd[1,4,lns[cur_wid][idx:length(lids[cur_wid])],
                        ] -= lvalue[cur_wid][idx:length(lids[cur_wid])]
                    end

                    idx = binarySearchFirstGreaterEqual(lids1[cur_wid], lft + 1)
                    if idx != -1
                        bs =lids1[cur_wid][idx:length(lids1[cur_wid])][app_rr[W[cur_wid][lids1[cur_wid][idx:length(
                                lids1[cur_wid],)]]]]
                        weight_upd[n24l, 4, W[cur_wid][bs]] += 1 ./ (rht .- bs)
                        change_id[W[cur_wid][bs]] .= 1
                    end
                    idx = binarySearchFirstLessEqual(rids[cur_wid], rht - 1)
                    weight_upd[1, 4, rns[cur_wid][1:idx]] -=
                        rvalue[cur_wid][1:idx]

                    idx = binarySearchFirstLessEqual(rids1[cur_wid], rht - 1)

                    bss=rids1[cur_wid][1:idx]
                    bs=bss[app_ll[W[cur_wid][bss]]]
                    weight_upd[n24r, 4, W[cur_wid][bs]] += 1 ./ (bs .- lft)
                    change_id[W[cur_wid][bs]] .= 1
                    idx = binarySearchFirstGreaterEqual(rids[cur_wid], rht + 1)
                    if idx != -1
                        ml=length(rns[cur_wid])
                        weight_upd[1, 4, rns[cur_wid][idx:ml]] -= rvalue[cur_wid][idx:ml]
                    end
                    if loc_l1[x][il] == loc_l1[y][jl] &&
                       loc_r1[x][ir] != loc_r1[y][jr]
                        tmp = loc_l2[x][il]
                        local_node = x
                        if W[cur_wid][lft] == x
                            tmp = loc_l2[y][jl]
                            local_node = y
                        end
                        v = 1 / (lenn - tmp)
                        weight_upd[1, 5-n24l, local_node] += v
                        tot[1, 5-n24l] += v
                    end
                    if loc_l1[x][il] != loc_l1[y][jl] &&
                       loc_r1[x][ir] == loc_r1[y][jr]
                        tmp = loc_r2[x][ir]
                        local_node = x
                        if W[cur_wid][rht] == x
                            tmp = loc_r2[y][jr]
                            local_node = y
                        end
                        v = 1 / (tmp - 1)
                        weight_upd[1, 5-n24r, local_node] += v
                        tot[1, 5-n24r] += v
                    end
                elseif W[cur_wid][lft] == x && W[cur_wid][rht] == y ||
                       W[cur_wid][lft] == y && W[cur_wid][rht] == x#1 x , 1 y
                    UN = unique(W[cur_wid][lft+1:rht-1])
                    v = 1 / (rht - lft)
                    weight_upd[2, 3, UN[xyinSam[UN]]] .-= v
                    totij += 1 / (p * (rht - lft))
                    lflag = !(W[cur_wid][rht] in W[cur_wid][2:We[cur_wid]])
                    rflag = !(W[cur_wid][lft] in W[cur_wid][We[cur_wid]+1:lenn-1])
                    if rflag
                        v = 1 / (lenn - lft)
                        tot[1, n24l] += v
                    end
                    if lflag
                        v = 1 / (rht - 1)
                        tot[1, n24r] += v
                    end
                    idx = binarySearchFirstLessEqual(lids[cur_wid], lft - 1)
                    weight_upd[1, 4, lns[cur_wid][1:idx]] -= lvalue[cur_wid][1:idx]
                    idx = binarySearchFirstGreaterEqual(lids[cur_wid], lft + 1)
                    if idx != -1
                        weight_upd[1,4,lns[cur_wid][idx:length(lids[cur_wid])],] -= lvalue[cur_wid][idx:length(lids[cur_wid])]
                    end


                    idx = binarySearchFirstGreaterEqual(lids1[cur_wid], lft + 1)
                    if idx != -1
                        bss=lids1[cur_wid][idx:length(lids1[cur_wid])]
                        bs=bss[app_rr[W[cur_wid][bss]]]
                        weight_upd[n24r, 4, W[cur_wid][bs]] += 1 ./ (rht .- bs)
                        change_id[W[cur_wid][bs]] .= 1
                    end
                    idx = binarySearchFirstGreaterEqual(rids[cur_wid], rht + 1)
                    if idx != -1
                        weight_upd[1,4,rns[cur_wid][idx:length(rids[cur_wid])],] -= rvalue[cur_wid][idx:length(rids[cur_wid])]
                    end
                    idx = binarySearchFirstLessEqual(rids[cur_wid], rht - 1)
                    weight_upd[1, 4, rns[cur_wid][1:idx]] -= rvalue[cur_wid][1:idx]

                    idx = binarySearchFirstLessEqual(rids1[cur_wid], rht - 1)
                    bss=rids1[cur_wid][1:idx]
                    bs=bss[app_ll[W[cur_wid][bss]]]
                    weight_upd[n24l, 4, W[cur_wid][bs]] += 1 ./ (bs .- lft)
                    change_id[W[cur_wid][bs]] .= 1
                end
            else
                if W[cur_wid][lft] == y && W[cur_wid][rht] != y#2 x, x==S, 1 y,y left
                    idx = binarySearchFirstLessEqual(lids[cur_wid], lft - 1)
                    weight_upd[1, 4, lns[cur_wid][1:idx]] -=
                        lvalue[cur_wid][1:idx]
                    UN = unique(W[cur_wid][lft+1:lenn-1])
                    v = 1 / (lenn - lft)
                    for t = 1:length(UN)
                        if xyinSam[UN[t]]
                            weight_upd[1, 3, UN[t]] -= v
                        end
                    end
                    weight_upd[1, 4, rns[cur_wid]] -= rvalue[cur_wid]
                    for b in rids1[cur_wid]
                        if app_ll[W[cur_wid][b]]
                            weight_upd[3, 4, W[cur_wid][b]] += 1 / (b - lft)
                            change_id[W[cur_wid][b]] = 1
                        end
                    end
                elseif W[cur_wid][lft] != y && W[cur_wid][rht] == y#2 x, x==S, 1 y,y right
                    idx = binarySearchFirstGreaterEqual(rids[cur_wid], rht + 1)
                    if idx != -1
                        weight_upd[1,4,rns[cur_wid][idx:length(rids[cur_wid])],
                        ] -= rvalue[cur_wid][idx:length(rids[cur_wid])]
                    end
                    UN = unique(W[cur_wid][2:rht-1])
                    v = 1 / (rht - 1)
                    for t = 1:length(UN)
                        if xyinSam[UN[t]]
                            weight_upd[1, 3, UN[t]] -= v
                        end
                    end
                    for b in lids1[cur_wid]
                        if app_rr[W[cur_wid][b]]
                            weight_upd[3, 4, W[cur_wid][b]] += 1 / (rht - b)
                            change_id[W[cur_wid][b]] = 1
                        end
                    end
                    weight_upd[1, 4, lns[cur_wid]] -= lvalue[cur_wid]
                elseif W[cur_wid][lft] == y && W[cur_wid][rht] == y#2 x, x==S, 2 y
                    weight_upd[1, 4, lns[cur_wid]] -= lvalue[cur_wid]

                    weight_upd[1, 4, rns[cur_wid]] -= rvalue[cur_wid]

                    idx = binarySearchFirstGreaterEqual(lids1[cur_wid], lft + 1)
                    if idx != -1
                        for b in lids1[cur_wid][idx:length(lids1[cur_wid])]
                            if app_rr[W[cur_wid][b]]
                                weight_upd[3, 4, W[cur_wid][b]] += 1 / (rht - b)
                                change_id[W[cur_wid][b]] = 1
                            end
                        end
                    end
                    idx = binarySearchFirstLessEqual(rids1[cur_wid], rht - 1)
                    for b in rids1[cur_wid][1:idx]
                        if app_ll[W[cur_wid][b]]
                            weight_upd[3, 4, W[cur_wid][b]] += 1 / (b - lft)
                            change_id[W[cur_wid][b]] = 1
                        end
                    end
                end
            end
            app_ll[W[cur_wid][lft:We[cur_wid]]] .= 1
            app_rr[W[cur_wid][rht:-1:(We[cur_wid]+1)]] .= 1
            il += loc_l1[x][il] == cur_wid
            ir += loc_r1[x][ir] == cur_wid
            jl += loc_l1[y][jl] == cur_wid
            jr += loc_r1[y][jr] == cur_wid
        end
        useless = 0
        change_id[x] = 1
        change_id[y] = 1
        change_id[S[1]] = 0
        changed_nodes = V[change_id]
        change_id[changed_nodes] .= 0
        ttt = zeros(Int, 3)
        lll = zeros(3, 3)
        ll = zeros(4, 4)
        weight_upd[:, :, changed_nodes] ./= p
        tot[1, 2] /= p
        tot[1, 3] /= p
        weight_upd[1, 1, changed_nodes] .= 0
        weight_upd[2, 2, changed_nodes] .= 0
        weight_upd[3, 3, changed_nodes] .= 0
        weight_upd[4, 4, changed_nodes] .= 0
        weight_upd[1, 2, changed_nodes] .+= cdct[x]
        weight_upd[1, 3, changed_nodes] .+= cdct[y]
        weight_upd[1, 4, changed_nodes] += cdct[changed_nodes]
        weight_upd[2, 3, changed_nodes] += weight_upd[3, 2, changed_nodes]
        weight_upd[2, 3, changed_nodes] .+= totij
        weight_upd[1, 2, changed_nodes] .-= tot[1, 2]
        weight_upd[1, 3, changed_nodes] .-= tot[1, 3]

        if x != S[1]
            for node in changed_nodes
                if node != x && node != y
                    weight_upd[2, 3, node] = weight_upd[2, 3, node]-1
                    t1=weight_upd[1, 2, node]
                    t2=weight_upd[1, 3, node]
                    t3=weight_upd[1, 4, node]
                    t4=weight_upd[2, 3, node]
                    t5=weight_upd[2, 4, node]
                    t6=weight_upd[3, 4, node]
                    ll[1, 2] = t1;oll[1,2] = t1;oll[2,1] = t1
                    ll[1, 3] = t2;oll[1,3] = t2;oll[3,1] = t2
                    ll[1, 4] = t3;oll[1,4] = t3;oll[4,1] = t3
                    ll[2, 3] = t4;oll[2,3] = t4;oll[3,2] = t4
                    ll[2, 4] = t5;oll[2,4] = t5;oll[4,2] = t5
                    ll[3, 4] = t6;oll[3,4] = t6;oll[4,3] = t6
                    if !isconnected(ll)
                        useless = 1
                        break
                    end
                    getER = getlittleER(ll,oll)
                    res_change[edge] += getER - iniRL[node]
                elseif node == x#2，4
                    weight_upd[1, 4, node] -= cdct[node]
                    weight_upd[1, 4, node] += weight_upd[1, 2, node]
                    weight_upd[3, 4, node] += weight_upd[2, 3, node]
                    weight_upd[3, 4, node] -= 1
                    t1=weight_upd[1, 3, node];
                    t2=weight_upd[1, 4, node];
                    t3=weight_upd[3, 4, node];
                    lll[1,2]=t1;olll[1,2]=t1;olll[2,1]=t1
                    lll[1,3]=t2;olll[1,3]=t2;olll[3,1]=t2
                    lll[2,3]=t3;olll[2,3]=t3;olll[3,2]=t3
                    if !isconnected(lll)
                        useless = 1
                        break
                    end
                    getER = getlittleER(lll,olll)
                    res_change[edge] += getER - iniRL[node]
                elseif node == y#3，4
                    weight_upd[1, 4, node] -= cdct[node]
                    weight_upd[1, 3, node] += weight_upd[1, 4, node]
                    weight_upd[2, 3, node] += weight_upd[2, 4, node]
                    weight_upd[2, 3, node] -= 1

                    for z = 1:3
                        for zz = z+1:3
                            t=weight_upd[z, zz, node]
                            lll[z, zz] = t
                            olll[z,zz] = t
                            olll[zz,z] = t
                        end
                    end

                    if !isconnected(lll)
                        useless = 1
                        break
                    end
                    getER = getlittleER(lll,olll)
                    res_change[edge] += getER - iniRL[node]
                end
            end
        else
            for node in changed_nodes
                if node != y
                    weight_upd[1, 3, node] += weight_upd[2, 3, node]
                    weight_upd[1, 4, node] += weight_upd[2, 4, node]
                    weight_upd[1, 3, node] -= 1
                    lll[1,2]= weight_upd[1, 3, node]
                    lll[1,3]= weight_upd[1, 4, node]
                    lll[2,3]= weight_upd[3, 4, node]
                    if !isconnected(lll)
                        useless = 1
                        break
                    end
                    getER = getlittleER(lll,olll)
                    res_change[edge] += getER - iniRL[node]
                elseif node == y
                    weight_upd[1, 4, node] -= cdct[node]
                    weight_upd[1, 3, node] += weight_upd[1, 4, node]
                    weight_upd[1, 3, node] += weight_upd[2, 3, node]
                    weight_upd[1, 3, node] += weight_upd[2, 4, node]
                    weight_upd[1, 3, node] -= 1
                    if weight_upd[1, 3, node] < 1e-6
                        useless = 1
                        break
                    end
                    res_change[edge] += (1 / weight_upd[1, 3, node]) - iniRL[node]
                end
            end
        end
        if useless == 1
            res_change[edge] = Inf
        end
        if !xsam
            xyinSam[x] = 0
        end
        if !ysam
            xyinSam[y] = 0
        end
    end
    xyinSam[sample_V] .= true
    res_change ./=sample_ratio
    res_change .+= app_inir
    for i = 1:G.m
        if res_change[i] < -1e8 || res_change[i] > 1e8
            res_change[i] = 0
        end
    end
    app_selc_edge = argmax(res_change)
    sclog = open("sclog.txt", "a")
    if mode=="bo"
        for i = 1:G.m
            if real_res[i] < -1e8 || real_res[i] > 1e8
                real_res[i] = 0
            end
        end
        real_selc_edge = argmax(real_res)
        logw(sclog, k, " select:",app_selc_edge, ",", real_selc_edge)
        logw(sclog, res_change[app_selc_edge], ",", real_res[real_selc_edge])
        del2[real_selc_edge] = 1
    else
        logw(sclog,k, " select:", app_selc_edge, "=",G.u[app_selc_edge]," ",G.v[app_selc_edge])
        logw(sclog, res_change[app_selc_edge])
    end
    close(sclog)
    del[app_selc_edge] = 1
    push!(deleted,app_selc_edge)
    xxx = G.u[app_selc_edge]
    yyy = G.v[app_selc_edge]
    deltt = 0
    nb1 = G.nbr[xxx]
    nbr_new = Int[]
    for i = 1:length(nb1)
        if deltt == 1 || nb1[i] != yyy
            push!(nbr_new, nb1[i])
        end
        if nb1[i] == yyy
            deltt = 1
        end
    end
    G.nbr[xxx] = nbr_new
    yyy = G.u[app_selc_edge]
    xxx = G.v[app_selc_edge]
    deltt = 0
    nb1 = G.nbr[xxx]
    nbr_new = Int[]
    for i = 1:length(nb1)
        if deltt == 1 || nb1[i] != yyy
            push!(nbr_new, nb1[i])
        end
        if nb1[i] == yyy
            deltt = 1
        end
    end
    G.nbr[xxx] = nbr_new
end
end
end
end
