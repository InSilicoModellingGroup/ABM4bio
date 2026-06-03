
# The ABM4bio project 

## What is this code about?

**ABM4bio** is a simple project that implements the agent-based modelling method to simulate systems related to biology, biomechanics and bioengineering.

## How to compile the code?

Before compiling the code of the **ABM4bio** project, you need to install **BioDynaMo** first. To do this, simply execute in your bash terminal the command `make biodynamo` which will take care of everything (i.e., updating dependencies, installing them in your computer, and finally compiling and installing **BioDynaMo** locally, i.e., in the *libs\\* directory of the **ABM4bio** project). Subsequently, you need to execute on the terminal the command `make fresh` that will compile the **ABM4bio** code and create the executable (available in the *build\\* directory of the project). To ensure everything is working in good order, you can execute on the terminal the command `make tests` taht will run one by one all simulation examples provided with this project.
