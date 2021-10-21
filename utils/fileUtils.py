import time

# not working need to fix
def verify_file_is_closed(path, cycles=10):
    cycle = 0
    while cycle < cycles:
        try:
            my_file = open(path, "r")
            if not my_file.buffer.raw.closefd:
                print("Please close file in path: {}".format(path))
                time.sleep(5)
                continue
            break
        except IOError:
            print("Please close file in path: {}".format(path))
            time.sleep(5)
