import Randomness_tests as tests
import miller_rabin as mr


def larger_closest_prime(x):
	i = 1
	while not mr.is_prime(x+i):
		i += 1
	return (x + i)

def smaller_closest_prime(x):
	if x < 3:
		return -100		# no smaller prime
	i = 1
	while not mr.is_prime(x-i):
		i += 1
	return (x - i)

def is_nearer_to_larger_closest_prime(x):
	larger = larger_closest_prime(x)
	smaller = smaller_closest_prime(x)
	distance_to_larger = larger - x
	distance_to_smaller = x - smaller
	if distance_to_larger < distance_to_smaller:	# return 1 if x is closer to the larger prime
		return 1
	elif distance_to_larger > distance_to_smaller:	# return -1 if x is closer to the smaller prime
		return -1
	else:						# return 0 if same distance
		return 0


# the generate_test_sequence function generate the sequence to be tested, 
# by storing the sequence in a list

# the return value is a list storing the test sequence

# can test on any sequence by modifying this function

# for the simplist case, can change the (num ** 2) in the parenthesis, 
# in order to generate some other simple sequence

# def generate_test_sequence(upper_bound):			# the test sequence is square numbers
# 	test_sequence = []
# 	num = 1
# 	index = 0
# 	while index < upper_bound:
# 		print('-Generating test sequence:{}%'.format(round(index*100 / upper_bound)), end='\r')
# 		test_sequence.append(num ** 2)				
# 		num += 1
# 		index += 1
# 	print()
# 	print('test sequence generated')
# 	print()
# 	return test_sequence

def generate_test_sequence(upper_bound):			# the test sequence is 2p where p are prime numbers
	test_sequence = []					# the corresponding bit sequence (generated with skip) appeared to be random
	num = 1
	index = 0
	while index < upper_bound:
		print('-Generating test sequence:{}%'.format(round(index*100 / upper_bound)), end='\r')
		num += 1
		if mr.is_prime(num):
			test_sequence.append(num * 2)
			index += 1
	print()
	print('test sequence generated')
	print()
	return test_sequence	

# generate the corresponding bit sequence from the chosen test sequence
# can change the 'skip' input (default 1) to avoid the correlation
def generate_bit_sequence(lower_bound, upper_bound, skip = 1):
	test_sequence = generate_test_sequence(upper_bound)
	bit_sequence = []
	for i in range(lower_bound, upper_bound, skip):
		print('-Generating bit sequence:{}%'.format(round(i*100 / (upper_bound - lower_bound))), end='\r')
		bit = is_nearer_to_larger_closest_prime(test_sequence[i])
		if bit != 0:	#ignore the same distance case
			bit_sequence.append(bit)

	for i in range(len(bit_sequence)):
		if bit_sequence[i] == -1:
			bit_sequence[i] = 0
	
	print()
	print('bit sequence generated')
	print()
	return bit_sequence



bit_sequence = generate_bit_sequence(100, 100000, 20)

# can choose any randomness tests combination from 'Randomness_tests.py'
# here are some examples
print('test result:')
tests.monobit(bit_sequence)
tests.serial(bit_sequence)
tests.poker(bit_sequence, 3)
tests.poker(bit_sequence, 4)
tests.runs(bit_sequence)
tests.autocor(bit_sequence, 4)



