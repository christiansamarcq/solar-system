# ========================================
# INTERACTIVE PYTHON TUTORIAL
# ========================================

def arithmetic_operations():
    print("\n=== ARITHMETIC OPERATIONS ===")
    print("Let's do some math!")

    while True:
        try:
            x = int(input("Enter first number: "))
            y = int(input("Enter second number: "))
            break
        except ValueError:
            print("Please enter whole numbers only (no decimals or letters)!")

    print(f"\nResults for {x} and {y}:")
    print(f"  Addition: {x} + {y} = {x + y}")
    print(f"  Subtraction: {x} - {y} = {x - y}")
    print(f"  Multiplication: {x} * {y} = {x * y}")
    print(f"  Division: {x} / {y} = {x / y}")
    print(f"  Integer Division: {x} // {y} = {x // y}")
    print(f"  Modulo (remainder): {x} % {y} = {x % y}")
    print(f"  Power: {x} ** {y} = {x ** y}")

    input("\nPress Enter to continue...")

def comparison_operations():
    print("\n=== COMPARISON OPERATIONS ===")
    print("These compare values and return True or False")

    while True:
        try:
            a = int(input("Enter first number: "))
            b = int(input("Enter second number: "))
            break
        except ValueError:
            print("Please enter whole numbers only!")

    print(f"\nComparing {a} and {b}:")
    print(f"  {a} == {b} (equal): {a == b}")
    print(f"  {a} != {b} (not equal): {a != b}")
    print(f"  {a} > {b} (greater): {a > b}")
    print(f"  {a} < {b} (less): {a < b}")
    print(f"  {a} >= {b} (greater or equal): {a >= b}")
    print(f"  {a} <= {b} (less or equal): {a <= b}")

    input("\nPress Enter to continue...")

def logical_operations():
    print("\n=== LOGICAL OPERATIONS ===")
    print("Combine conditions with AND, OR, NOT")

    print("\nExample: Can you drive?")
    while True:
        try:
            age = int(input("How old are you? "))
            break
        except ValueError:
            print("Please enter a whole number for age!")
    has_license = input("Do you have a license? (yes/no): ").lower() == "yes"

    is_adult = age >= 18
    can_drive = is_adult and has_license

    print(f"\n  Age >= 18: {is_adult}")
    print(f"  Has license: {has_license}")
    print(f"  Can drive (both True): {can_drive}")

    input("\nPress Enter to continue...")

def for_loops():
    print("\n=== FOR LOOPS ===")
    print("Repeat code a specific number of times\n")

    print("1. Simple counting")
    while True:
        try:
            n = int(input("Count from 1 to what number? "))
            break
        except ValueError:
            print("Please enter a whole number!")
    for i in range(1, n + 1):
        print(i, end=" ")
    print()

    print("\n2. Multiplication table")
    while True:
        try:
            num = int(input("Which multiplication table (1-12)? "))
            break
        except ValueError:
            print("Please enter a whole number!")
    for i in range(1, 11):
        print(f"{num} x {i} = {num * i}")

    print("\n3. Loop through your items")
    items = input("Enter items separated by commas: ").split(",")
    for item in items:
        print(f"  - {item.strip()}")

    input("\nPress Enter to continue...")

def while_loops():
    print("\n=== WHILE LOOPS ===")
    print("Repeat while a condition is True\n")

    print("1. Countdown")
    while True:
        try:
            start = int(input("Count down from what number? "))
            break
        except ValueError:
            print("Please enter a whole number!")
    while start > 0:
        print(start, end=" ")
        start -= 1
    print("Blast off!")

    print("\n2. Number guessing game")
    secret = 7
    guess = 0
    attempts = 0

    print("I'm thinking of a number between 1 and 10...")
    while guess != secret:
        while True:
            try:
                guess = int(input("Your guess: "))
                break
            except ValueError:
                print("Please enter a whole number!")
        attempts += 1
        if guess < secret:
            print("Too low!")
        elif guess > secret:
            print("Too high!")
        else:
            print(f"Correct! You got it in {attempts} attempts!")

    input("\nPress Enter to continue...")

def loop_control():
    print("\n=== LOOP CONTROL (break & continue) ===")

    print("1. BREAK - Exit loop early")
    print("Finding numbers divisible by 7:")
    while True:
        try:
            limit = int(input("Search up to what number? "))
            break
        except ValueError:
            print("Please enter a whole number!")
    for i in range(1, limit + 1):
        if i % 7 == 0:
            print(f"Found: {i}")
            if input("Keep searching? (yes/no): ").lower() != "yes":
                break

    print("\n2. CONTINUE - Skip to next iteration")
    print("Printing only even numbers:")
    while True:
        try:
            max_num = int(input("Print even numbers from 1 to: "))
            break
        except ValueError:
            print("Please enter a whole number!")
    for i in range(1, max_num + 1):
        if i % 2 != 0:  # If odd
            continue  # Skip odd numbers
        print(i, end=" ")
    print()

    input("\nPress Enter to continue...")

def functions_tutorial():
    print("\n=== FUNCTIONS ===")
    print("Functions are reusable blocks of code that do specific tasks\n")

    # Example 1: Simple function
    print("1. BASIC FUNCTION")
    print("Let's create a greeting function!\n")

    name = input("What's your name? ")

    # Define the function
    def greet(person_name):
        print(f"Hello, {person_name}! Welcome to Python!")
        print(f"Nice to meet you, {person_name}!")

    print(f"\nCalling greet('{name}'):")
    greet(name)

    # Example 2: Function with return value
    print("\n2. FUNCTIONS THAT RETURN VALUES")
    print("Functions can calculate and return results\n")

    while True:
        try:
            num1 = int(input("Enter first number: "))
            num2 = int(input("Enter second number: "))
            break
        except ValueError:
            print("Please enter whole numbers!")

    # Define calculator functions
    def add(a, b):
        return a + b

    def multiply(a, b):
        return a * b

    def power(a, b):
        return a ** b

    print(f"\nResults:")
    print(f"  add({num1}, {num2}) = {add(num1, num2)}")
    print(f"  multiply({num1}, {num2}) = {multiply(num1, num2)}")
    print(f"  power({num1}, {num2}) = {power(num1, num2)}")

    # Example 3: Function that checks conditions
    print("\n3. FUNCTIONS WITH DECISION MAKING")
    print("Functions can use if/else to make decisions\n")

    while True:
        try:
            age = int(input("Enter an age: "))
            break
        except ValueError:
            print("Please enter a whole number!")

    def check_age_category(age):
        if age < 13:
            return "Child"
        elif age < 20:
            return "Teenager"
        elif age < 65:
            return "Adult"
        else:
            return "Senior"

    category = check_age_category(age)
    print(f"Age {age} is in category: {category}")

    # Example 4: Build your own function
    print("\n4. CREATE YOUR OWN FUNCTION")
    print("Let's build a custom calculator function!\n")

    while True:
        try:
            x = int(input("Enter a number: "))
            break
        except ValueError:
            print("Please enter a whole number!")

    print("\nWhat operation?")
    print("  1. Square it (x * x)")
    print("  2. Double it (x * 2)")
    print("  3. Half it (x / 2)")
    operation = input("Choose 1, 2, or 3: ")

    def custom_calculator(number, op):
        if op == "1":
            result = number * number
            return f"{number} squared = {result}"
        elif op == "2":
            result = number * 2
            return f"{number} doubled = {result}"
        elif op == "3":
            result = number / 2
            return f"{number} halved = {result}"
        else:
            return "Invalid operation"

    print(f"\nResult: {custom_calculator(x, operation)}")

    # Example 5: Function with default parameters
    print("\n5. FUNCTIONS WITH DEFAULT VALUES")
    print("Parameters can have default values if not provided\n")

    def create_greeting(name, greeting="Hello", punctuation="!"):
        return f"{greeting}, {name}{punctuation}"

    your_name = input("Enter your name: ")

    print("\nDifferent ways to call the same function:")
    print(f"  Default: {create_greeting(your_name)}")
    print(f"  Custom greeting: {create_greeting(your_name, 'Hi')}")
    print(f"  Custom all: {create_greeting(your_name, 'Hey', '!!!')}")

    # Summary
    print("\n" + "="*40)
    print("KEY CONCEPTS:")
    print("="*40)
    print("• def - keyword to create a function")
    print("• Parameters - inputs to the function")
    print("• return - send a value back")
    print("• Calling - using the function by name()")
    print("• Default values - optional parameters")
    print("="*40)

    input("\nPress Enter to continue...")

def main_menu():
    while True:
        print("\n" + "="*40)
        print("PYTHON TUTORIAL - Choose a topic:")
        print("="*40)
        print("1. Arithmetic Operations (+, -, *, /, etc.)")
        print("2. Comparison Operations (==, >, <, etc.)")
        print("3. Logical Operations (and, or, not)")
        print("4. For Loops")
        print("5. While Loops")
        print("6. Loop Control (break & continue)")
        print("7. Functions (def, return, parameters)")
        print("8. Exit")
        print("="*40)

        choice = input("Enter your choice (1-8): ")

        if choice == "1":
            arithmetic_operations()
        elif choice == "2":
            comparison_operations()
        elif choice == "3":
            logical_operations()
        elif choice == "4":
            for_loops()
        elif choice == "5":
            while_loops()
        elif choice == "6":
            loop_control()
        elif choice == "7":
            functions_tutorial()
        elif choice == "8":
            print("\nHappy coding!")
            break
        else:
            print("Invalid choice. Please try again.")

# Start the interactive tutorial
if __name__ == "__main__":
    print("Welcome to the Interactive Python Tutorial!")
    main_menu()
