// C program for the above approach
#include <stdio.h>
typedef enum {
  QR_NAIVE = 0, // mostly previous PPRZ
  QR = 1,
  CHOL = 2,
  CG = 3,
  } activeSetAlgoChoice;

typedef activeSetAlgoChoice (*activeSetVariant)(int y[2], const float);
 
activeSetAlgoChoice fun1(int* y, int x)
{
    return *y + 10 - x;
}
activeSetAlgoChoice fun2(int* y, int x)
{
    return *y + 30 - x;
}
 
// Function that return type ptr
activeSetVariant fun(int choice)
{
    if (choice == 0)
        return &fun1;
    else
        return &fun2;
}
 
// Driver Code
int main()
{
    int a = 10;
 
    printf("%d\n", fun(0)(&a, 1));
    printf("%d\n", fun(1)(&a, 2));
 
    return 0;
}