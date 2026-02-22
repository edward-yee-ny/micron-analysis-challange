import math, csv
from collections import defaultdict
from itertools import product as iproduct

MINS_PW = 7 * 24 * 60
UTIL = {'A':0.78,'B':0.76,'C':0.80,'D':0.80,'E':0.76,'F':0.80,
        'A+':0.84,'B+':0.81,'C+':0.86,'D+':0.88,'E+':0.84,'F+':0.90}
SPACE = {'A':6.78,'B':3.96,'C':5.82,'D':5.61,'E':4.65,'F':3.68,
         'A+':6.93,'B+':3.72,'C+':5.75,'D+':5.74,'E+':4.80,'F+':3.57}
FAB_LIM = {'Fab1':1500,'Fab2':1300,'Fab3':700}
ALL_WS = ['A','B','C','D','E','F','A+','B+','C+','D+','E+','F+']

nodes = {
    'Node1': [(1,'D',14,'D+',12),(2,'F',25,'F+',21),(3,'F',27,'F+',23),
              (4,'A',20,'A+',16),(5,'F',12,'F+',9),(6,'D',27,'D+',21),
              (7,'D',17,'D+',13),(8,'A',18,'A+',16),(9,'A',16,'A+',13),
              (10,'D',14,'D+',11),(11,'F',18,'F+',16)],
    'Node2': [(1,'F',19,'F+',16),(2,'B',20,'B+',18),(3,'E',10,'E+',7),
              (4,'B',25,'B+',19),(5,'B',15,'B+',11),(6,'F',16,'F+',14),
              (7,'F',17,'F+',15),(8,'B',22,'B+',16),(9,'E',7,'E+',6),
              (10,'E',9,'E+',7),(11,'E',20,'E+',19),(12,'F',21,'F+',18),
              (13,'E',12,'E+',9),(14,'E',15,'E+',12),(15,'E',13,'E+',10)],
    'Node3': [(1,'C',21,'C+',20),(2,'D',9,'D+',7),(3,'E',24,'E+',23),
              (4,'E',15,'E+',11),(5,'F',16,'F+',14),(6,'D',12,'D+',11),
              (7,'C',24,'C+',21),(8,'C',19,'C+',13),(9,'D',15,'D+',13),
              (10,'D',24,'D+',20),(11,'E',17,'E+',15),(12,'E',18,'E+',13),
              (13,'F',20,'F+',18),(14,'C',12,'C+',11),(15,'D',11,'D+',10),
              (16,'C',25,'C+',20),(17,'F',14,'F+',13)],
}
loading = {'Node1':[12000,10000,8500,7500,6000,5000,4000,2000],
           'Node2':[5000,5200,5400,5600,6000,6500,7000,7500],
           'Node3':[3000,4500,7000,8000,9000,11000,13000,16000]}
quarters = ["Q1'26","Q2'26","Q3'26","Q4'26","Q1'27","Q2'27","Q3'27","Q4'27"]

N3_C_STEPS = {s[0] for s in nodes['Node3'] if s[1]=='C'}

def compute_tools_needed(n1, n2, n3, f1_n1, f1_n2,
                          n3_f1_c, n3_f2_c,   # Node3 C-steps in Fab1, Fab2
                          n3_f1_nc, n3_f2_nc, n3_f3_nc):
    """
    Compute exact min tools needed per fab per WS type.
    Uses mintech RPTs. Returns {fab: {ws: count}}.
    ALL fabs can have any WS type.
    """
    mins = {fab: defaultdict(float) for fab in ['Fab1','Fab2','Fab3']}
    f2_n1 = n1 - f1_n1
    f2_n2 = n2 - f1_n2
    n3_f3_c = n3 - n3_f1_c - n3_f2_c

    for s, ws, rpt, wst, rptt in nodes['Node1']:
        mins['Fab1'][ws] += f1_n1 * rpt
        mins['Fab2'][ws] += f2_n1 * rpt
    for s, ws, rpt, wst, rptt in nodes['Node2']:
        mins['Fab1'][ws] += f1_n2 * rpt
        mins['Fab2'][ws] += f2_n2 * rpt
    for s, ws, rpt, wst, rptt in nodes['Node3']:
        if s in N3_C_STEPS:
            mins['Fab1'][ws] += n3_f1_c * rpt
            mins['Fab2'][ws] += n3_f2_c * rpt
            mins['Fab3'][ws] += n3_f3_c * rpt
        else:
            mins['Fab1'][ws] += n3_f1_nc * rpt
            mins['Fab2'][ws] += n3_f2_nc * rpt
            mins['Fab3'][ws] += n3_f3_nc * rpt

    tools = {fab:{} for fab in ['Fab1','Fab2','Fab3']}
    for fab in ['Fab1','Fab2','Fab3']:
        for ws, m in mins[fab].items():
            if m > 0:
                tools[fab][ws] = math.ceil(m / (MINS_PW * UTIL[ws]))
    return tools

def get_space(tools):
    return {fab: sum(cnt * SPACE[ws] for ws, cnt in tools[fab].items() if ws in SPACE)
            for fab in ['Fab1','Fab2','Fab3']}

def space_ok(tools):
    sp = get_space(tools)
    return all(sp[f] <= FAB_LIM[f] for f in ['Fab1','Fab2','Fab3']), sp

all_flows = []
all_tool_plans = []

print("="*65)
print("FINAL FEASIBLE SOLVER")
print("Tools are per-quarter; C-steps can go to any fab")
print("="*65)

for qi, q in enumerate(quarters):
    n1 = loading['Node1'][qi]
    n2 = loading['Node2'][qi]
    n3 = loading['Node3'][qi]

    found = None

    # Search: r1 = Node1 Fab1 frac, r2 = Node2 Fab1 frac
    # n3_nc_f1 = Node3 NC Fab1 frac (rest in Fab2, none in Fab3 for simplicity)
    # n3_c_f3 = Node3 C Fab3 frac (rest split Fab1/Fab2)
    best_space_sum = float('inf')

    ratios = [i/10 for i in range(11)]
    for r1, r2, nc_f1_r, c_f3_r in iproduct(ratios, ratios, ratios, ratios):
        f1_n1 = round(n1 * r1)
        f1_n2 = round(n2 * r2)
        n3_f1_nc = round(n3 * nc_f1_r)
        n3_f2_nc = n3 - n3_f1_nc
        n3_f3_c = round(n3 * c_f3_r)
        n3_rem_c = n3 - n3_f3_c
        n3_f1_c = n3_rem_c // 2
        n3_f2_c = n3_rem_c - n3_f1_c

        tools = compute_tools_needed(n1, n2, n3, f1_n1, f1_n2,
                                     n3_f1_c, n3_f2_c,
                                     n3_f1_nc, n3_f2_nc, 0)
        ok, sp = space_ok(tools)
        if ok:
            total_sp = sum(sp.values())
            if total_sp < best_space_sum:
                best_space_sum = total_sp
                found = (r1, r2, nc_f1_r, c_f3_r,
                         f1_n1, n1-f1_n1, f1_n2, n2-f1_n2,
                         n3_f1_nc, n3_f2_nc, 0,
                         n3_f1_c, n3_f2_c, n3_f3_c,
                         tools, sp)

    if found is None:
        print(f"{q}: NO FEASIBLE SOLUTION FOUND - space constraints impossible!")
        # Fallback: minimize Fab3 C load at the cost of violating space
        f1_n1 = n1//2; f1_n2 = n2//2
        n3_f1_nc = n3//2; n3_f2_nc = n3 - n3_f1_nc
        n3_f3_c = 0; n3_f1_c = n3//2; n3_f2_c = n3 - n3_f1_c
        tools = compute_tools_needed(n1,n2,n3,f1_n1,f1_n2,n3_f1_c,n3_f2_c,n3_f1_nc,n3_f2_nc,0)
        _, sp = space_ok(tools)
        found = (0.5,0.5,0.5,0,f1_n1,n1-f1_n1,f1_n2,n2-f1_n2,
                 n3_f1_nc,n3_f2_nc,0,n3_f1_c,n3_f2_c,0,tools,sp)

    (r1,r2,nc_f1_r,c_f3_r,f1n1,f2n1,f1n2,f2n2,
     n3f1nc,n3f2nc,n3f3nc,n3f1c,n3f2c,n3f3c,
     tools,sp) = found

    ok_str = "OK" if all(sp[f]<=FAB_LIM[f] for f in sp) else "SPACE VIOLATION!"
    print(f"\n{q}: {ok_str}")
    print(f"  N1={n1} (F1={f1n1}/F2={f2n1}) | N2={n2} (F1={f1n2}/F2={f2n2})")
    print(f"  N3={n3}: C-steps F1={n3f1c}/F2={n3f2c}/F3={n3f3c} | NC F1={n3f1nc}/F2={n3f2nc}")
    print(f"  Space: F1={sp['Fab1']:.0f}/1500  F2={sp['Fab2']:.0f}/1300  F3={sp['Fab3']:.0f}/700")
    for fab in ['Fab1','Fab2','Fab3']:
        items = {ws:cnt for ws,cnt in tools[fab].items() if cnt>0}
        if items: print(f"  {fab} tools: {items}")

    # Build flow
    flow = {'Node1':{},'Node2':{},'Node3':{}}
    for s,ws,rpt,wst,rptt in nodes['Node1']:
        flow['Node1'][s] = {'Fab1':f1n1,'Fab2':f2n1,'Fab3':0}
    for s,ws,rpt,wst,rptt in nodes['Node2']:
        flow['Node2'][s] = {'Fab1':f1n2,'Fab2':f2n2,'Fab3':0}
    for s,ws,rpt,wst,rptt in nodes['Node3']:
        if s in N3_C_STEPS:
            flow['Node3'][s] = {'Fab1':n3f1c,'Fab2':n3f2c,'Fab3':n3f3c}
        else:
            flow['Node3'][s] = {'Fab1':n3f1nc,'Fab2':n3f2nc,'Fab3':0}

    all_flows.append({'q':q,'flow':flow})
    all_tool_plans.append({'q':q,'tools':tools})

# Write flow CSV
with open('/home/claude/flow_ok.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['Quarter','Node','Step','Fab','Loading (to fill)'])
    for r in all_flows:
        for node in ['Node1','Node2','Node3']:
            nn = node.replace('Node','')
            for s,*_ in nodes[node]:
                d = r['flow'][node][s]
                w.writerow([r['q'],nn,s,1,d['Fab1']])
                w.writerow([r['q'],nn,s,2,d['Fab2']])
                w.writerow([r['q'],nn,s,3,d['Fab3']])

# Write tools CSV - EXACT tools needed per quarter
with open('/home/claude/tools_ok.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['Quarter','WS Name','Fab 1 - WS Count','Fab 2 - WS Count','Fab 3 - WS Count'])
    for r in all_tool_plans:
        for ws in ALL_WS:
            w.writerow([r['q'], ws,
                        r['tools']['Fab1'].get(ws,0),
                        r['tools']['Fab2'].get(ws,0),
                        r['tools']['Fab3'].get(ws,0)])

print("\n\nCSVs written: flow_ok.csv and tools_ok.csv")

python3 -c "
import math
# Q2'27-Q4'27 still infeasible. Let me compute exactly why.
# Q2'27: n3=11000. Node3 C-steps need certain space, Node1+Node2 need certain space.
# With C-steps split across Fab1+Fab2 (not Fab3), what's the minimum space?

MINS_PW = 7*24*60
UTIL = {'A':0.78,'B':0.76,'C':0.80,'D':0.80,'E':0.76,'F':0.80}
SPACE = {'A':6.78,'B':3.96,'C':5.82,'D':5.61,'E':4.65,'F':3.68}

# Minimum machines needed in total for all nodes combined in Q2'27:
# Node1=5000, Node2=6500, Node3=11000
# Node1 uses: A,D,F
# Node2 uses: B,E,F
# Node3 uses: C,D,E,F (non-C) and C (C-steps)

# If we put everything in one fab with infinite space, what's the total?
nodes1_rpt = {'D':14+27+17+14,'F':25+27+12+18,'A':20+18+16}  # sum RPTs per WS
nodes2_rpt = {'B':20+25+15+22,'E':10+7+9+20+12+15+13,'F':19+16+17+21}
nodes3c_rpt = {'C':21+24+19+12+25}  # C-step RPTs
nodes3nc_rpt = {'D':9+12+15+24+11,'E':24+15+17+18,'F':16+20+14}

n1=5000; n2=6500; n3=11000

total_mins = {}
for ws,rpt in nodes1_rpt.items():
    total_mins[ws] = total_mins.get(ws,0) + n1*rpt
for ws,rpt in nodes2_rpt.items():
    total_mins[ws] = total_mins.get(ws,0) + n2*rpt
for ws,rpt in nodes3c_rpt.items():
    total_mins[ws] = total_mins.get(ws,0) + n3*rpt
for ws,rpt in nodes3nc_rpt.items():
    total_mins[ws] = total_mins.get(ws,0) + n3*rpt

total_tools = {ws: math.ceil(m/(MINS_PW*UTIL[ws])) for ws,m in total_mins.items()}
total_space = sum(t*SPACE[ws] for ws,t in total_tools.items())

print('Q2 27 total tools if single fab:')
for ws,t in sorted(total_tools.items()): print(f'  {ws}: {t} tools = {t*SPACE[ws]:.0f}m2')
print(f'Total space needed: {total_space:.0f}m2')
print()
print(f'Available space across 3 fabs: 1500+1300+700=3500m2')
print(f'Total needed: {total_space:.0f}m2')
print(f'Fits in total space: {total_space <= 3500}')

python3 -c "
import math
# Even across ALL 3 fabs total (3500m2), the tools needed for Q2'27 = 3606m2.
# This is 106m2 over. This is a FUNDAMENTAL infeasibility in the problem itself.
# 
# The ONLY escape: use TOR tools which have SMALLER footprint per tool.
# TOR tools are more efficient per wafer (lower RPT) so fewer tools needed.
# AND some TOR tools are smaller (e.g., F+=3.57 vs F=3.68, D+=5.74 vs D=5.61)
# 
# Let's check with ALL TOR tools:
MINS_PW = 7*24*60
UTIL = {'A+':0.84,'B+':0.81,'C+':0.86,'D+':0.88,'E+':0.84,'F+':0.90}
SPACE = {'A+':6.93,'B+':3.72,'C+':5.75,'D+':5.74,'E+':4.80,'F+':3.57}

# TOR RPTs
nodes1_rpt = {'D+':12+21+13+11,'F+':21+23+9+16,'A+':16+16+13}
nodes2_rpt = {'B+':18+19+11+16,'E+':7+6+7+19+9+12+10,'F+':16+14+15+18}
nodes3c_rpt = {'C+':20+21+13+11+20}
nodes3nc_rpt = {'D+':7+11+13+20+10,'E+':23+11+15+13,'F+':14+18+13}

n1=5000; n2=6500; n3=11000
total_mins = {}
for ws,rpt in nodes1_rpt.items():
    total_mins[ws] = total_mins.get(ws,0) + n1*rpt
for ws,rpt in nodes2_rpt.items():
    total_mins[ws] = total_mins.get(ws,0) + n2*rpt
for ws,rpt in nodes3c_rpt.items():
    total_mins[ws] = total_mins.get(ws,0) + n3*rpt
for ws,rpt in nodes3nc_rpt.items():
    total_mins[ws] = total_mins.get(ws,0) + n3*rpt

total_tools = {ws: math.ceil(m/(MINS_PW*UTIL[ws])) for ws,m in total_mins.items()}
total_space = sum(t*SPACE[ws] for ws,t in total_tools.items())

print('Q2 27 ALL TOR tools:')
for ws,t in sorted(total_tools.items()): print(f'  {ws}: {t} tools = {t*SPACE[ws]:.0f}m2')
print(f'Total: {total_space:.0f}m2 / 3500m2 available -> fits: {total_space <= 3500}')
print()
# Check for Q3'27 and Q4'27
for n1,n2,n3,q in [(4000,7000,13000,\"Q3'27\"),(2000,7500,16000,\"Q4'27\")]:
    total_mins2 = {}
    for ws,rpt in nodes1_rpt.items(): total_mins2[ws]=total_mins2.get(ws,0)+n1*rpt
    for ws,rpt in nodes2_rpt.items(): total_mins2[ws]=total_mins2.get(ws,0)+n2*rpt
    for ws,rpt in nodes3c_rpt.items(): total_mins2[ws]=total_mins2.get(ws,0)+n3*rpt
    for ws,rpt in nodes3nc_rpt.items(): total_mins2[ws]=total_mins2.get(ws,0)+n3*rpt
    t2 = {ws: math.ceil(m/(MINS_PW*UTIL[ws])) for ws,m in total_mins2.items()}
    s2 = sum(t*SPACE[ws] for ws,t in t2.items())
    print(f'{q} ALL TOR: {s2:.0f}m2 / 3500m2 -> fits: {s2<=3500}')