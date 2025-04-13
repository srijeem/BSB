import sys

def get_mapping(map_file):
    # Open the file.
    f = open(map_file, "r")

    # Result is a list of dictionaries.
    mapping_list = []

    # Skip the header on the first line.
    header = f.readline()

    for line in f:
        #TODO: PUT YOUR CODE HERE
        lines = line.strip('\n').split('\t')
        #no of dictionaries created
        d = len(lines) - 1
        # accessing the Ensemble ID 
        e_id = lines[0]

        
        #checking whether the no of dicts in the list is equals to d or not
        if len(mapping_list) != d:
            # if thats not the case, then append no of dicts to mapping_list which will be equals to d
            i = 0
            
            while i <= (d):
                mapping_list.append({})
                i+=1
            
        
        # loop through the non-Ensembl IDs
        # j=1
        # while j<dict+1:
        for j in range(1, d+1):
        
            # skip empty elements in lines
            if lines[j] == '':
                continue
            # check if non-Ensembl and Ensembl ID are already present in the dict of mapping_list
            elif (lines[j], e_id) not in mapping_list[j - 1]:
                # add non-Ensembl ID and Ensembl ID as key,value to the dict in mapping_list
                mapping_list[j - 1][lines[j]] = e_id
            # j+=1


    # Remember to close the file after we're done.
    f.close()

    return mapping_list


def get_go_terms(mapping_list, go_file):
    # Open the file.
    f = open(go_file, "r")

    # This will be the dictionary that this function returns.
    # Entries will have as a key an Ensembl ID and the value will
    # be a set of GO terms.
    go_dict = dict()

    for line in f:
        # TODO : PUT YOUR CODE HERE
        if line[0] == '!':
            continue
        else:
            
            e_id = ''
            
            lines = line.split('\t')
            # get the IDs of the protiens from lines
            p_id = lines[1]

            # iterating through the dictionaries inside the mapping_list
            for k in mapping_list:
            
                # check whether the (protein ID) p_id is a key or not
                if p_id in k:
                # then look at the Ensembl ID in dict and add it to the (Ensembl ID) e_id
                    e_id = k[p_id]
                # otherwise continue(skip) if  the (protein_id) p_id not present in dict
                else:
                    continue
            # if the(protein ID) p_id if not found in dictionaries inside the mapping_list, then skip
            if e_id == '':
                continue

           # # set go_term equal to 4th item in lines and remove 'GO:'
            go_term = lines[4].replace('GO:', '')

            # check whether the (Ensembl ID) e_id is a key in the go_dict dictionary
            if e_id not in go_dict:
            #     # if not then, add it as key and set go_term as its value
                go_dict[e_id] = {go_term}
            else:
            #     # otherwise if the e_id is a key, add go_term to set
                go_dict[e_id].add(go_term)
            # go_dict[e_id]={go_term} if  e_id not in go_dict else go_dict[e_id].add(go_term)



    # Remember to close the file after we're done.
    f.close()

    return go_dict


def compute_score(alignment_file, go_one_dict, go_two_dict):
    # Open the file.
    f = open(alignment_file, "r")

    # Keep track of the number of proteins we can't map to GO terms
    # and the score.
    unmappable_one = 0
    unmappable_two = 0
    score = 0.0

    for line in f:
        #TODO : PUT YOUR CODE HERE
        for line in f:
        ## create list out of e_id list split by tabs
            lines = line.split()
        #if first e_id is not in go_one_dict, then add 1 to unmappable_one
            if lines[0] not in go_one_dict:
                unmappable_one += 1
            #if second e_id not in go_two_dict, then add 1 to unmappable_two
            if lines[1] not in go_two_dict:
                unmappable_two += 1
            # if both the e_ids are in it,then, calculate the jaccard index and add the value to the variable score 
            if lines[0] in go_one_dict:
                if lines[1] in go_two_dict:
                    score=score + len(go_one_dict[lines[0]] & go_two_dict[lines[1]]) / len(
                    go_one_dict[lines[0]] | go_two_dict[lines[1]])
                #     score=score+ len(go_one_dict[lines[0]]).intersection(go_two_dict[lines[1]]) / len(
                # go_one_dict[lines[0]]).union(go_two_dict[lines[1]])


    # Remember to close the file after we're done.
    f.close()

    # Return the statistics and the score back so the main code
    # can print it out.
    return unmappable_one, unmappable_two, score


def main():
    # TODO: PUT YOUR CODE HERE

    #calling the mapping function
    map1 = get_mapping(sys.argv[4])
    map2 = get_mapping(sys.argv[5])

    # calling the get_go_terms function 
    go1 = get_go_terms(map1, sys.argv[2])
    go2 = get_go_terms(map2, sys.argv[3])
    
    #Writting Intermediate Output Files
    with open(sys.path[0]+'/map1.txt','a') as outmap1:
        outmap1.write(str(map1)+'\n')
    with open(sys.path[0]+'/map2.txt','a') as outmap2:
        outmap2.write(str(map2)+'\n')
    with open(sys.path[0]+'/go1.txt','a') as outgo1:
        outgo1.write(str(go1)+'\n')
    with open(sys.path[0]+'/go2.txt','a') as outgo2:
        outgo2.write(str(go2)+'\n')

    # check if passed no of arguments are correct or not
    if len(sys.argv) != 6:
        sys.exit('ERROR! The number of arguments passed is INCORRECT! You need to enter 5 arguments.')

    print('\nScoring Network Alignments Script Executed\n')
    #writting intermediate files)
    print('Below Intermediate Output Files are generated:\n')
    print('map1.txt, map2.txt, go1.txt, go2.txt\n')
    # call compute_score function and print return values
    print('Computed Score:\t',compute_score(sys.argv[1], go1, go2))


if __name__ == '__main__':
    main()