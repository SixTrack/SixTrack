#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestPlugin.h"

int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}

