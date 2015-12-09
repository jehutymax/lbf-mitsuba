CS613 - Machine Learning

To run this implementation, which contains:
 - my contributions
 - reference filter implementation from original work
 - pbrt renderer

 you'll need:
  - Visual Studio 2013
  - CUDA 6.5 (version is important - CUDA 7.5 contains its own half type which has nothing to do with Ilmbase's half class, and then there's a conflict since namespaces are not being used, but should)


  Open the solution in Visual Studio 2013:
   lbf-pbrt\src\solution.vs2013\pbrt.sln

   And hit Build Solution.
   All 10 projects should build successfully.


   Then, open the Visual Studio 2013 x64 Command Prompt, and navigate to lbf-pbrt\pbrt-scenes

   Then, in the command-line:
   .\..\..\pbrt.exe {any pbrt file in the folder} --spp {number of samples}

   All images should be output to the same working directory.