import glob
import random


def read_names():
    """
    Read the names in ../names
    :return: list of names
    """
    names = []
    for name_file in glob.glob("../names/*"):
        with open(name_file, 'r') as f:
            first_line = f.readline()
            first_name = first_line.split()[0]
            names.append(first_name)
    return names


greetings = ["Hello", "Hi", "Howdy", "Hola", "Ciao", "G'day", "What's up", "Howzgarn"]


def greet_random(names):
    global greetings
    print(f"{random.choice(greetings)}, {random.choice(names)}!")


def greet_all(names):
    global greetings
    for name in names:
        print(f"{random.choice(greetings)}, {name}!")


def greet_each_other(names):
    global greetings
    name1, name2 = random.sample(names, 2)

    print(f"{random.choice(greetings)}, my name is {name1}, what's your name?")
    print(f"{random.choice(greetings)}, my name is {name2}, How are you enjoying this course?")
    print("I am loving it, but I wish the conference was in-person and we could all be there...")


def main():
    names = read_names()
    greet_random(names)
    greet_all(names)
    greet_each_other(names)


if __name__ == '__main__':
    main()
