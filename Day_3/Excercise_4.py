# 🚨 Don't change the code below 👇
print("Welcome to Python Pizza Deliveries!")
size = input("What size pizza do you want? S, M, or L ")
add_pepperoni = input("Do you want pepperoni? Y or N ")
extra_cheese = input("Do you want extra cheese? Y or N ")
# 🚨 Don't change the code above 👆

#Write your code below this line 👇
small = 15
medium = 20
large = 25
spep = 2
mlpep = 3
xcheese = 1
sum = 0
if size == "S":
    sum += small
elif size == "M":
    sum += medium
elif size == "L":
    sum += large

if add_pepperoni == "Y":
    if size == "S":
        sum += spep
    else:
        sum += mlpep

if extra_cheese == "Y":
    sum += xcheese

print(f"Your final bill is: ${sum}.")