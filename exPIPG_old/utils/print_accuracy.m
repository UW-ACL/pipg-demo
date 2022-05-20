function rel_err = print_accuracy(sol1,sol2)
    rel_err = 100*norm(sol1.Z-sol2.Z)/norm(sol2.Z);
    fprintf("\n\n%s rel. accuracy wrt %s : %.2e %%\n",sol1.name,sol2.name,rel_err);
end