import re
import os

files = [
    "site/components/4_results.html",
    "site/components/5_discussion.html",
    "site/components/6_conclusion.html"
]

replacements = [
    (r'\(step_11\)', r''),
    (r'\(step_12; ', r'('),
    (r'\(step_12\)', r''),
    (r'\(step_13\)', r''),
    (r'\(step_13, ', r'('),
    (r'\(step_14\)', r''),
    (r'\(step_19\)', r''),
    (r'Yes \(step_19\)', r'Yes'),
    (r'\(step_05 ', r'('),
    (r'\(step_15\)', r''),
    (r'\(step_17\)', r''),
    (r'\(step_06, ', r'('),
    (r'Step_16 reports', r'A correlated significance synthesis reports'),
    (r'Step_16 then reports', r'A correlated significance synthesis is then reported.'),
    (r'Step_16 correlated', r'Correlated'),
    (r'Step_15 external-informed', r'External-informed'),
    (r'Step_17 directional-odds', r'Directional-odds'),
    (r'Step_17 adds', r'An additional analysis adds'),
    (r'step_15 applies', r'an external-informed test applies'),
    (r'step_11', r'the model-dependence analysis'),
    (r'step_12', r'the microlensing-nuisance Monte Carlo'),
    (r'step_13', r'the hierarchical Bayesian comparison'),
    (r'Step_14 established', r'Initial integration established'),
    (r'Step_19 then provides', r'A subsequent standalone ingestion then provides'),
    (r'First-pass ingestion \(step_14\) pulled', r'First-pass ingestion pulled'),
    (r'standalone expansion step \(step_19\) then', r'standalone expansion step then'),
    (r'step_14', r'the initial integration'),
    (r'step_19', r'the subsequent standalone ingestion'),
    (r'step_05', r'a prior step'),
    (r'step_15', r'an external-informed test'),
    (r'step_17', r'a complementary analysis'),
    (r'step_06', r'a prior sensitivity scan'),
    (r'step_16', r'a synthesis step'),
    (r'step_04', r''),
    (r'step_07', r''),
    (r'step_08', r''),
    (r'step_09', r''),
    (r'step_10', r''),
]

for filename in files:
    filepath = os.path.join("/Users/matthewsmawfield/www/TEP-LENS", filename)
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        continue
        
    with open(filepath, "r") as f:
        content = f.read()
    
    for old, new in replacements:
        content = re.sub(old, new, content)
        
    # Clean up empty parens and double spaces
    content = content.replace(" ()", "")
    content = content.replace("( ", "(")
    content = content.replace("  ", " ")
    
    with open(filepath, "w") as f:
        f.write(content)

print("Done")
