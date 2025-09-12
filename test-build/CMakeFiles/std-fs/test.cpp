
#include <filesystem>
#include <iostream>

int main(int argc, char* argv[]){
  std::cout << std::filesystem::temp_directory_path();
  return 0;
}