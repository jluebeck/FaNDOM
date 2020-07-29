def parse_cmap(cmapf, keep_length=False):
    cmaps = {}
    # contigCovs = {}
    with open(cmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                if fD["CMapId"] not in cmaps:
                    cmaps[fD["CMapId"]] = {}
                    # contigCovs[fD["CMapId"]] = {}

                # this is not a good way to parse label channel means color channel
                if fD["LabelChannel"] == "1":
                    cmaps[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Position"])
                    # contigCovs[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Coverage"])

                elif fD["LabelChannel"] == "0" and keep_length:
                    cmaps[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Position"])
    return cmaps


def parse_bnx(bnxF, keep_length=False):
    moleculeD = {}
    with open(bnxF) as infile:
        for line in infile:
            if not line.startswith('#'):
                fields = line.rstrip().rsplit("\t")
                if line.startswith('0'):
                    currKey = fields[1]
                elif line.startswith('1'):
                    if keep_length:
                        moleculeD[currKey] = [float(x) for x in fields[1:]]
                    else:
                        moleculeD[currKey] = [float(x) for x in fields[1:-1]]

    return moleculeD


def parse_xmap(xmapf):
    detailFields = ["XmapEntryID", "QryContigID", "RefContigID", "Orientation", "Confidence", "QryLen", "RefLen",
                    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Alignment"]

    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                alnstring = ")" + fD["Alignment"] + "("
                # xmapAln[fD["XmapEntryID"]] = alnstring.rsplit(")(")[1:-1]

                xmapPair[fD["XmapEntryID"]] = {x: fD[x] for x in detailFields}

    return xmapPair
contig_id_map = {1: 1, 12: 2, 16: 3, 17: 4, 18: 5, 19: 6, 20: 7, 21: 8, 22: 9, 2: 10, 3: 11, 4: 12, 5: 13, 6: 14,
                 7: 15,
                 8: 16, 9: 17, 10: 18, 11: 19, 13: 20, 14: 21, 15: 22, 24: 23, 25: 24}

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", help="query directory", required=True)
parser.add_argument("-r", "--ref", help="reference directory", required=True)
parser.add_argument("-a", "--align", help="FaNDOM directory", required=True)
# parser.add_argument("-o", "--output", help="Output directory", required=True)
args = parser.parse_args()

ref = parse_cmap(args.ref)

xmap = parse_xmap(args.align)
if args.query.endswith('.cmap'):
    mol = parse_cmap(args.query)
elif args.query.endswith('.bnx'):
    mol = parse_bnx(args.query)
    
for k in xmap.keys():
    # print(k)
    a = xmap[k]
    q_id = a['QryContigID']
    ref_id = a['RefContigID']
    alignment = a['Alignment']
    alignment = alignment.split(')')[:-1]
    for i in range(0, len(alignment) - 1):
        ref_1 = alignment[i].split(',')[0][1:]
        q_1 = alignment[i].split(',')[1]
        ref_2 = alignment[i + 1].split(',')[0][1:]
        q_2 = alignment[i + 1].split(',')[1]
        # print(ref_id, ' a ', q_id)
        ref_dist = abs(ref[ref_id][int(ref_1)] - ref[ref_id][int(ref_2)])
        q_dist = abs(mol[q_id][int(q_1)] - mol[q_id][int(q_2)])
        if q_dist - ref_dist > 3000:
            print('insertion' + '\t' + str(q_id) + '\t' + str(contig_id_map[int(ref_id)]) + '\t' + str(mol[q_id][int(q_1)]) + '\t' + str(mol[q_id][
                int(q_2)]) + '\t' + str(ref[ref_id][
                      int(ref_1)]) + '\t' + str(ref[ref_id][
                      int(ref_2)]) + '\t' + str(q_1) + '\t' + q_2 + '\t' + ref_1 + '\t' + ref_2+'\t'+ str(abs(q_dist - ref_dist)) )
        elif ref_dist - q_dist > 3000:
            print('deletion' + '\t' + q_id + '\t' + str(contig_id_map[int(ref_id)]) + '\t' + str(mol[q_id][int(q_1)]) + '\t' + str(mol[q_id][
                int(q_2)]) + '\t' + str(ref[ref_id][
                      int(ref_1)]) + '\t' + str(ref[ref_id][
                      int(ref_2)]) + '\t' + q_1 + '\t' + q_2 + '\t' + ref_1 + '\t' + ref_2 +'\t'+str(abs(q_dist - ref_dist)))
