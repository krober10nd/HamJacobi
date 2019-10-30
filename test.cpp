#include<iostream> 
#include<algorithm> 
#include<vector> 

int main() {

  std::vector<int> a {1,2,3,4,1,6};
  
  auto iter = a.begin();
  while ((iter = std::find(iter, a.end(), 1)) != a.end())
  {
      std::cout<<std::distance(a.begin(),iter)<<std::endl;
      iter++;
  }

  return 0; 
}

