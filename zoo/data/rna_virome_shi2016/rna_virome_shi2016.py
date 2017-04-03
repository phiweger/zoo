acc = ['KX' + str(i) for i in range(882764, 884872 + 1)]
selection = acc[:100]
with open('rna_virome_shi2016.txt', 'w+') as file:
    for i in selection:
        file.write(''.join([i, '\n']))
# len 2109
