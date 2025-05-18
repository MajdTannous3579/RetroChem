import json
import os

path = input("Enter 'json' file: ").strip()
if not os.path.isfile(path):
    print("Invalid path!")
    exit(1)

ret = []
desc = []
with open(path) as file:
    reactions = json.load(file)
    for reaction in reactions:
        ret.append(reactions[reaction]['retro_smarts'].replace('\\', '-').replace('/', '-'))
        desc.append(reactions[reaction]['description'])

with open("datamol.db", "w") as file:
    file.write('[\n')
    if len(ret) != 0:
        file.write(f'  [\n    "{ret[0]}",\n    {{"description" : "{desc[0]}"}}\n  ]')
    for i in range(1, len(ret)):
        file.write(f',\n  [\n    "{ret[i]}",\n    {{"description" : "{desc[i]}"}}\n  ]')
    file.write('\n]\n')