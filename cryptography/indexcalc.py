import random
import time
# fast exponentiation
def power(a, x, n):
  ret = 1
  if x == 0:
    return ret % n
  if (x % 2 == 1):
    ret = a
  return (ret * (power(a, x // 2, n)**2)) % n

def gcd(a, b):
    while(b):
        a, b = b, a%b
    return a

def hasInverse(n, mod):
    if(gcd(n, mod)!=1):
        return False
    return True

def modInverse(a, m) :
    """
    Function for finding modInverse
    """
    m0 = m 
    y = 0
    x = 1
    if (m == 1) : 
        return 0
    while (a > 1) :

        q = a // m 
        t = m 

        m = a % m 
        a = t 
        t = y 
  
        # Update x and y 
        y = x - q * y 
        x = t 
  
    # Make x positive 
    if (x < 0) : 
        x = x + m0 
    return x
def modInverse2(a, m):
    """
    Function for finding modInverse assuming gcd(a,m) = 1 and m is prime
    """
    return power(a, m-2, m)    
def torref(M, order):
    """
    Function for finding the Reduced Row Echelon Form of matrix M (mod order)
    """
    if not M:
        return
    lead = 0
    rowCount = len(M)
    columnCount = len(M[0])-1
    for r in range(rowCount):
        #rref found
        if lead >= columnCount:
            return
        #Begin with leftmost nonzero column
        i = r
        
        #Get the first row with nonzero value in the pivot column
        while(M[i][lead]==0):
            i += 1
            #No nonzero row found
            if(i==rowCount):
                i = r
                #Start row reduction from next column
                lead += 1
                #rref found
                if(columnCount==lead):
                    return
        #Swap the row to the pivot
        M[i],M[r] = M[r],M[i]
        
        #Use row addition operations to create zeros in all positions below the pivot
        pivot = M[r][lead]
        
        #Check if the pivot has an inverse to do "divide" with
        #Possible failure if all values in the column do not have an inverse
        if(hasInverse(pivot, order)):
            inv = modInverse(pivot, order)
            
            #Reduce the pivot to 1 by "dividing"
            M[r] = [ mrx*inv%order for mrx in M[r]]
            
            #Set the rest of the values in the column to 0
            for i in range(rowCount):
                if (i!= r):
                    pivot = M[i][lead]
                    M[i] = [ (iv - pivot*rv)%order for rv,iv in zip(M[r],M[i])]
            lead += 1
        #Get the next row with nonzero value

def removeZero(M):
    """
    Function for getting matrix M excluding rows with zeroes in the augment col
    """
    rowCount = len(M)
    col = len(M[0])-1
    newM =[]
    for r in range(rowCount):
        if(M[r][col]!=0):
            newM.append(M[r])
    return newM
# assume n is prime
# since b ^ (n - 1) = 1(mod n)
# a): b ^ s = 1(mod n)
# b): b ^ s*(2^ri) = -1(mod n) for some ri in [0, r)
# if b) is not true, then b ^ s*(2^ri) = - 1(mod n) 
# we can keep take square root of b ^ s*(2^ri) until ri = 0,
# at that time we have b ^ s = 1(mod n)

def miller_rabin_test(b, r, s, n):
  # 1: test if b ^ s = 1(mod n) or -1(mod n)
  b_s = power(b, s, n)
  if b_s == 1 or b_s == n-1:
    return True

  # 2: test if b ^ s*(2^ri) = -1(mod n) for some ri in [0, r)
  for i in range(r - 1):
    b_s = (b_s * b_s) % n
    if b_s == (n - 1):
      return True

  return False

def miller_rabin(n, k):
  # easy case
  if n == 2 or n == 3 or n == 5:
    return True
  if n == 4:
    return False

  # step 1: find r,s where s * 2^r = n - 1
  r = 0
  s = n - 1;
  while s > 0 and (s % 2 == 0):
    r = r + 1
    s = s // 2

  # step 2: randomly chose base b, do k times miller_rabin primality test
  for i in range(k):
    if not miller_rabin_test(random.randint(2, n - 2),r, s, n):
      return False
  return True

def factor_base_hit(n, factor_base):
  relation = [0] * len(factor_base)

  # factor n to our factor base
  for i in range(0, len(factor_base)):
    if n == 1: break

    while n % factor_base[i] == 0:
      relation[i] = relation[i] + 1
      n = n // factor_base[i]

  if n == 1:
    return relation
  else:
    return []

def figure_out_index(relations, factor_base, relation, offset):
  x = -offset
  for i in range(0, len(factor_base)):
    x = x + (relations[i][-1] * relation[i])
  return x

def getRelations(g, p, factor_base,relations, slack, smoothness):
  while True:
    k = random.randint(2, p-2)
    # test if we already have enough system relations
    if len(relations) >= (slack + smoothness):
      break

    # initialize potential relation
    g_k = power(g, k, p)
    relation = factor_base_hit(g_k, factor_base)

    # if g_k completely factor to our base primes record the relation
    if relation:
      relation.append(k)
      relations.append(relation)
def isRREF(relations):
    a = 1
    for i in range(len(relations)):
        a=a*relations[i][i]
    if(a==1):
        return True
    return False

def isQR(n, p):
    """
    Function for checking if n is a quadratic residue of odd p
    Return:
      1 if n is a quadratic residue
      0 if n = 0(mod p)
      -1 if n is not a quadratic residue
    """
    return power(n, (p-1)//2, p)

# solve g^x = h(mod p)
def index_calculus(g, h, p, smoothness, slack, order):
  pstart = time.time()
  # step 0: find out all primes that are smaller than
  factor_base = [2]
  i = 3
  while True:
    if len(factor_base) >= smoothness:
      break
    #If we use 3, check isQR(i, p) here
    if miller_rabin(i, 100):
      factor_base.append(i)
    i = i + 2
  ftime = time.time()
  print("Time to get factor bases: "+str(ftime-pstart))
  print ('successfully initialize' , smoothness, 'factor bases')

  # step 1: find system relations for g^k
  
  relations = []
  while((len(relations)!=smoothness) or not isRREF(relations)):
      start = time.time()
      getRelations(g, p, factor_base,relations, slack, smoothness)
      time1 = time.time()
      print("Time to get relations: "+str(time1-start))
      torref(relations, order)
      relations = removeZero(relations)
      time2 = time.time()
      print("Time to get rref: "+str(time2-time1))
      print ('successfully get',  len(relations), ' relations')
      for rw in relations:
          print (" ".join( (str(rv) for rv in rw)))
  
  #step 2: try to factorize h*g^s into our factor base
  while True:
    k = random.randint(2, p-2)
    g_k_h = (power(g, k ,p) * h) % p
    factor_base_hits = factor_base_hit(g_k_h, factor_base)
    if factor_base_hits:
      x = figure_out_index(relations, factor_base, factor_base_hits, k)%(order)
      time3 = time.time()
      print("Time to solve: "+str(time3-time2))
      print ('Solved -> ', x)
      break;

def main():
 #  index_calculus(5, 20, 503, 10, 30, 502)
 #  405
  index_calculus(2, 92327518017225, 247457076132467, 500, 500, 247457076132466)
  # 208891284998759
#   index_calculus(2, 31, 83, 10, 15)

if __name__ == "__main__":
  main()
