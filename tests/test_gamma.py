import json
with open('results/outputs/step_03_tep_closure.json') as f:
    d = json.load(f)
print(d['tep_closure_loops']['S1_S4_SX'])
