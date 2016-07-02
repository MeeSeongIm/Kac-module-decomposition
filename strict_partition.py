from itertools import chain, combinations
import collections


def powerset(seq):
    if len(seq) <= 1:
        yield seq
        yield []
    else:
        for item in powerset(seq[1:]):
            yield [seq[0]]+item
            yield item

# Notation: trans_fun(d) = trans_fun(n) for the Lie superalgebra p(n)
# Change the number 2 to any n (you also have to type in n not necessarily distinct integers for the program to work). 
# scroll to the bottom of the script

def pwr(d,r):
    out = []
    for n in range(2**d):
        p = [int(r)*int(i) for i in list(bin(n)[2:])]
        while len(p) < d:
            p = [0] + p
        out.append(p)
    return out
    print(out)


 

def trans_fun(d):
    yes = set(['yes', 'y', 'ye', 's'])
    no = set(['no', 'n'])
    if d<=0:
        print ("error; d must be a positive integer. Run this program again.")
    else:
        print ("Let d be a positive integer. Wait for the next command.")
        arr = [0]*d # initial array of length d has only 0's in its entries 
        for i in range(len(arr)):
        	arr[i] = int(input("Partition (integer) entry (from left to right) and then hit ENTER (each time): ")) # don't hit ENTER without typing in a number.
        print ("Your partition is ", arr,".")
        print(" ")
        bigKac = arr # make a copy of arr as bigKac

        transFunctor0 = [arr for x in range(d)] # initial partition multiple times before any application of the translation functor 
        applyTransFunctor = input("Apply translation functor S (y or n)?")
        applyTransFunctor = applyTransFunctor.lower()
        elem = [[0 for x in range(d)] for x in range(d)] # will be an array with e_i's
        negElem = [[0 for x in range(d)] for x in range(d)] # will be an array with -e_i's
        for i in range(d):
            elem[i][i] = elem[i][i]+1
            negElem[i][i] = negElem[i][i]-1
        basicElem = elem # an array of distinct elem basis vectors
        negBasicElem = negElem
        lengthTwoBasisElem = elem
        lengthTwoNegBasisElem = negElem
        for i in range(d):
            lengthTwoBasisElem[i][i] = elem[i][i]+1
            lengthTwoNegBasisElem[i][i] = negElem[i][i]-1
        bigKacTObigKac = []
        

        PowerSetOfD = [[] for _ in range(1)]

        NumbersUpToD = range(1,int(d)+1,1)
        SetOfD = set(NumbersUpToD)
        Numbers = [i for i in range(1)]

        for z in chain.from_iterable(combinations(SetOfD,r) for r in range(len(SetOfD)+1)):
        	PowerSetOfD[Numbers[0]].append(z) 
        PowerSetOfD = sum(PowerSetOfD, [])
        #print(PowerSetOfD) # prints the power set of positive integers 1 to d: {1,...,d}

        powerSetCopiesbigKac = [arr for x in range(2**d)] # power set of copies of big Kac 
        lengTwoBasisElem = [[0 for i in range(d)] for i in range(2**d)]
        
        for i in range(d):
        	lengTwoBasisElem[i][i] = lengTwoBasisElem[i][i] + 2
        
        #r = [x for x in powerset(lengTwoBasisElem)]

        
        r = pwr(d,2)
        
        #for i in range(2**d):
        #	for j in range(d):
        #		r[i][j] = int(0 if value is None else value)
        	
        lengTwoBasisElem = r
        
        #print(lengTwoBasisElem)

 
 
 
        if applyTransFunctor in no:
            print("Your partition is", arr,". The code has been written to stop here.")
        elif applyTransFunctor in yes:
            transFunctor1 = [[0 for x in range(d)] for x in range(d)]
            transFunctorNeg1 = [[0 for x in range(d)] for x in range(d)]
            print(" ")
            howManyTimes_TransFunctor = int(input("How many times? Positive integer only (this response is the last one needed to completely run the rest of this program): "))
            #if howManyTimes_TransFunctor in set():
            #	print("Apply tranlation functor 0 times requested.")
            #   quit

            for j in range(howManyTimes_TransFunctor):
                for i in range((2**j)*(d**(j+1))):
                    transFunctor1[i] = [x+y for x,y in zip(transFunctor0[i],elem[i])]
                    transFunctorNeg1[i] = [x+y for x,y in zip(transFunctor0[i],negElem[i])]
                transitionFunctor0 = transFunctor1 + transFunctorNeg1
                print(" ")
                print("Applying the translation functor", j+1,"time(s), we obtain the following set of weights: ",transitionFunctor0)
                print("There are", len(transitionFunctor0), "(not necessarily distinct) weights when applying the translation functor", j+1,"time(s).")
                #print("checking previous step's elem:", elem) # for checking purposes 
                #print("checking previous step's negElem:", negElem) # for checking purposes 
                prevElem = elem
                prevNegElem = negElem
                elem = [[basicElem[m] for x in range((2**(j+1))*(d**(j+1)))] for m in range(len(basicElem))] # make elem to have enough copies of elem vectors
                negElem =[[negBasicElem[m] for x in range((2**(j+1))*(d**(j+1)))] for m in range(len(negBasicElem))]
                elem = sum(elem, []) # analogous to flattening
                negElem = sum(negElem, [])
                bigKacTObigKac = transitionFunctor0
                transitionFunctor0 = [transitionFunctor0 for u in range(d)] # weight vector d times (no need to fix)
                transitionFunctor0 = sum(transitionFunctor0, [])
                transFunctor0 = transitionFunctor0 # renaming to transFunctor0 
                transFunctor1 = [[0 for x in range(d)] for x in range(len(elem))]
                transFunctorNeg1 = [[0 for x in range(d)] for x in range(len(negElem))]
                #print(transFunctor0) # print d copies of all the weights
                #print(len(transFunctor0))



                #Decompose the big Kac module in terms of small Kac modules.
                smallKac = [[0 for x in range(d)] for x in range(2**d)]
                for i in range(2**d):
                  smallKac[i] = [x+y for x, y in zip(powerSetCopiesbigKac[i] ,lengTwoBasisElem[i])]
            print(" ")
            print("The big Kac module decomposes into the following small Kac modules whose partitions are:",smallKac)
            print(" ")
            
            bigKacTObigKac = list(map(lambda x: tuple(x),bigKacTObigKac))
            res1 = dict()
            
            for i in bigKacTObigKac:
            	if i not in res1:
            		res1[i] = 1
            	else:
            		res1[i] +=1
            
            print("The big Kac module with partition", bigKac, "tensored with V numerous times has decomposition into the following big Kac modules whose partitions are given, together with their multiplicities:",res1)
            
            smallKac = list(map(lambda x: tuple(x), smallKac))
            res2 = dict()

            print(" ")	            
            for i in smallKac:
            	if i not in res2:
            		res2[i] = 1
            	else:
            		res2[i] += 1
            		
            print("The big Kac module with partition", bigKac, "decomposes into the following small Kac modules whose partitions are given, together with their multiplicities:",res2)
            	
            commonWeightsBigKac = {x:res1[x] for x in res1 if x in res2}
            commonWeightsSmallKac = {x:res2[x] for x in res1 if x in res2}
            print(" ")
            print("Partitions and their multiplicities from the decomposition as big Kac modules that also appear as partitions as small Kac:", commonWeightsBigKac) 
            print(" ")
            print("Partitions and their multiplicities from the decomposition as small Kac modules that also appear as partitions as big Kac:",commonWeightsSmallKac)
            
            multiplicity = {x:commonWeightsBigKac[x]*commonWeightsSmallKac[x] for x in commonWeightsBigKac}
            print(" ")
            print("Conclusion: thus the partitions appearing, together with their multiplicities, in the endomorphism algebra invariant under p(d) are:" ,multiplicity)



        else:
        	while (applyTransFunctor not in no) or (applyTransFunctor not in yes):
        	    print("You've answered no. End of the code.")
        	    break

# Notation: trans_fun(d) = trans_fun(n) for the Lie superalgebra p(n)
# Change the number 2 to any n (you also have to type in n not necessarily distinct integers for the program to work). 

trans_fun(2) 








# trans_fun()
# d = int(input("Enter a number for d:"))
 

# repeat for the small Kac module, combine the partitions that they have in common. 


# want to type in n ordered numbers
# apply the translation functor S = sum S_i  on big Kac module d (or 2d) times (add e_i and subtract e_i one coordinate at a time and print each output, say A, this is a list), call this S^d (big Kac module of weight lambda)
# produce the tuples and save as txt file
# apply the translation functor S ??? times on a small Kac module (add 2 e_i one coordinate at a time and print each output, say B, this is the final list of weights)
# compare A and B: fix items in A and go through items in B, if there are any in common, then produce a copy onto list C (count multiplicities)
# make a dictionary D: {({tuple: multiplicity},{tuple: multiplicity},{}}
#
# command line: first we enter the weight vector lambda
# type: all, prints dictionary
# type specific weight: produce multiplicity
# type dimension: produce the sum of the multiplicities
# the dimension is the dim of endom of p(n)-invariant subalgebra









