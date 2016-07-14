

manager = DBManager("savedData_0/SingleParamVariation.json")


paramter = ["alpha, beta"]
definedRange = (range[0,10, range[1,2]])

for paramter in range(definedRange):
    manager.record_Data(alpha=alpha, beta=0, tau=0, gamma=0, time=t, timestep=dt,outputData =
                       simulation.run(alpha, beta, tau, gamma, time, timestep))


# after all sims have been run
manager.save_Data()