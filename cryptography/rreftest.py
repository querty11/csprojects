
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
                if i != r:
                    pivot = M[i][lead]
                    M[i] = [ (iv - pivot*rv)%order for rv,iv in zip(M[r],M[i])]
            lead += 1

        #Get the next row with nonzero value
        else:
            continue
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
def Legendre(n, p):
    """
    Function for checking if n is a quadratic residue of odd prime p
    Return:
      1 if n is a quadratic residue
      -1 if n is not a quadratic residue
    """
    return power(n, (p-1)//2, p)

if __name__ == '__main__':
    #We should choose suitable values for testing for g, p, and order
    g = 5
    p = 503
    order = p-1
    m = [[  0 , 1 , 0 , 1 , 0 , 0 , 1 , 0 , 0 , 0, 471], 
[  6 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 208], 
[  0 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 0 , 0, 353], 
[  3 , 1 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0, 497], 
[  0 , 3 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0,  90], 
[  2 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0,  59], 
[  0 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0, 224], 
[  0 , 3 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0, 195], 
[  0 , 2 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 314], 
[  3 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0, 190], 
[  1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 1, 428], 
[  1 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0, 439], 
[  0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 1 , 0, 141], 
[  1 , 0 , 1 , 0 , 0 , 1 , 0 , 0 , 0 , 0, 327], 
[  0 , 0 , 2 , 0 , 0 , 0 , 1 , 0 , 0 , 0, 231], 
[  0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 0, 149], 
[  0 , 1 , 1 , 0 , 1 , 0 , 0 , 0 , 0 , 0, 199], 
[  1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0, 114], 
[  2 , 0 , 0 , 2 , 0 , 0 , 0 , 0 , 0 , 0,  74], 
[  0 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0, 224], 
[  0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 1, 349], 
[  1 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0, 326], 
[  3 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0, 190], 
[  7 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0,  64], 
[  6 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 209], 
[  1 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 0 , 0,  53], 
[  1 , 0 , 0 , 1 , 0 , 0 , 0 , 1 , 0 , 0,  23], 
[  0 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 157], 
[  3 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 260], 
[  0 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 0, 466], 
[  0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 0 , 0, 474], 
[  2 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0, 490], 
[  1 , 0 , 0 , 0 , 0 , 1 , 0 , 1 , 0 , 0,  61], 
[  0 , 2 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 314], 
[  0 , 2 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 0, 399], 
[  0 , 1 , 0 , 1 , 0 , 0 , 0 , 1 , 0 , 0, 479], 
[  0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0, 124], 
[  2 , 2 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 215], 
[  3 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 260], 
[  1 , 0 , 0 , 1 , 0 , 0 , 1 , 0 , 0 , 0,  15]]
    torref(m, order)
    m = removeZero(m)
    for rw in m:
        print (" ".join( (str(rv) for rv in rw)))
    
