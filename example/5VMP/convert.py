import parmed as pmd
amber = pmd.load_file('complex_wat.prmtop',xyz='complex_wat.inpcrd')
amber.save('complex_wat.gro')
amber.save('complex_wat.top')
quit()

