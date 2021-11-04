  auto start = std::chrono::steady_clock::now();
  run_ea(argc, argv, ea, fit_t());//fit_t() is optionnal
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  // basic file operations

//still BAD as it also times the writing into the files

  std::ofstream myfile;
  myfile.open ("timing.dat");
  myfile << "\n Time elapsed (s): " << elapsed_seconds.count();
  myfile << std::endl;
  myfile.close();
