import os, json

# import json
# with open('Hjets/333_H1jB_card.json','r') as f:
#     j = json.load(f)
#     print j['output']

# with open('Hjets/333_H1jB_card.json','r+') as f:
#     j = json.load(f)
#     # print j['input'][0]['files'][0]
#     input = j['input'][0]
#     files = input['files']
#     print files[0]
#     input['original_files'] = files[:]
#     for i in range(len(files)):
#         files[i] = '$TMPDIR/'+os.path.basename(files[i])
#     # input['files'] = [ '$TMPDIR/'+os.path.basename(x) for x in files ]
#     print files[0]
#     f.seek(0)
#     json.dump(j,f)
#     f.truncate()

with open('Hjets/333_H1jB_card.json','r+') as f:
    j = json.load(f)
    input = j['input'][0]
    files = input['files']
    print files[0]
    input['original_files'] = files[:]
    for i in range(len(files)):
        files[i] = '$TMPDIR/'+os.path.basename(files[i])
    print files[0]
    f.seek(0)
    json.dump(j,f)
    f.truncate()

