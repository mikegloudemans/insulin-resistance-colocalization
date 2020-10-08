import requests


with open("data/curated_gene_sets/single_coloc_genes.txt") as f:
    with open("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/mouse_knockouts/drug_targets.txt", "w") as w:
        for line in f:
            gene = line.strip()
            r = requests.get('https://platform-api.opentargets.io/v3/platform/public/search?q={0}'.format(gene))
            output = r.json()

            found_match = False
            for i in range(len(output['data'])):
                if 'approved_symbol' not in output['data'][i]['data']:
                    continue

                if gene != output['data'][i]['data']['approved_symbol']:
                    print "No match with", gene, output['data'][i]['data']['approved_symbol']
                    continue
                found_match = True  
                break

            if found_match:

                print gene
                print output['data'][i]['data']['drugs']
                w.write("{0}\t{1}\n".format(gene, output['data'][i]['data']['drugs']))
            else:
                w.write("{0}\t{1}\n".format(gene, "no-match"))
                
