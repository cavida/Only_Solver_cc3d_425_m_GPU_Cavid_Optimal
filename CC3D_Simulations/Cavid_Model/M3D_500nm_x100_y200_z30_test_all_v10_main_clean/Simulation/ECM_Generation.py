from cc3d.core.PySteppables import *
import random
random.seed(123)


class ECM_Generation(SteppableBasePy):

    def __init__(self, frequency):

        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.direction_dict = {1: (+1, 0, 0),  # those starting with 0X like 1,2,3,... are lower directions
                          2: (+1, +1, 0),
                          3: (0, +1, 0),
                          4: (-1, +1, 0),
                          5: (-1, 0, 0),
                          6: (-1, -1, 0),
                          7: (0, -1, 0),
                          8: (+1, -1, 0),

                          11: (+1, 0, -1),  # those starting with 1X like 11,12,13,... are flat plane direction
                          12: (+1, +1, -1),
                          13: (0, +1, -1),
                          14: (-1, +1, -1),
                          15: (-1, 0, -1),
                          16: (-1, -1, -1),
                          17: (0, -1, -1),
                          18: (+1, -1, -1),

                          21: (+1, 0, 1),  # those starting with 2X like 21,22,23,... are uppper directions
                          22: (+1, +1, 1),
                          23: (0, +1, 1),
                          24: (-1, +1, 1),
                          25: (-1, 0, 1),
                          26: (-1, -1, 1),
                          27: (0, -1, 1),
                          28: (+1, -1, 1),

                          40: (0, 0, 1),  # upward direction
                          30: (0, 0, -1)}  # downward direction

        x_limit = self.dim.x
        self.y_limit = self.dim.y - 30
        z_limit = self.dim.z
        self.all_directions = tuple(self.direction_dict.keys())
        basic_direction_weights = (1, 1, 1, 1,  # Down direction weigths
                                   1, 1, 1, 1,  # Flat plane direction weight
                                   1, 1, 1, 1,  # Up direction weights
                                   1)  # Pole to pole direction weight (z)

        self.w = {}
        for i in self.all_directions:
            if i < 10:
                if i <= 4:
                    self.w[i] = basic_direction_weights[i - 1]
                else:
                    self.w[i] = basic_direction_weights[i - 1 - 4]

            elif 10 < i < 20:
                if i <= 14:
                    self.w[i] = basic_direction_weights[i - 1 - 10 + 4]
                else:
                    self.w[i] = basic_direction_weights[i - 1 - 10 + 4 - 4]

            elif 20 < i < 30:
                if i <= 24:
                    self.w[i] = basic_direction_weights[i - 1 - 20 + 8]
                else:
                    self.w[i] = basic_direction_weights[i - 1 - 20 + 8 - 4]

            else:
                self.w[i] = basic_direction_weights[12]
        print('&&&&&&&',self.dim.x)
        print('********************', self.w)

        num_Fibers = (1 + int(15 * self.dim.x * self.y_limit /(500 * 500))) * self.dim.z
        for j in range(num_Fibers):  # Insert the number for different percentage of collagen fibers (17)
            x = random.randint(0, self.dim.x)
            y = random.randint(0, self.y_limit)  # CHANGED THIS to 19 FOR ENDO
            z = random.randint(0, self.dim.z)
            # direction = random.randint(1, 26) #Random directions
            init_direction = random.choices(self.all_directions, self.w)[0]
            # direction = 1
            Fiber_Length = 20  # pixels fiber length

            #This method call is used for generating each fiber; it accepts the initial location, initial direction, and fiber length
            self.ECM_Gen(x,y,z,init_direction,Fiber_Length,self.ECM)

        print("ECM Generation is Finished!")

    def step(self, mcs):
        Fiber_Length = 20  #fiber length
        for cell in self.cell_list_by_type(self.AVEC):
            print('Cell',len(self.cell_list_by_type(self.AVEC)))
            pixel_list = self.get_cell_boundary_pixel_list(cell)
            r = random.choice(list(pixel_list))  # Choosing a random pixel at the boundary of a cell
            init_direction = random.choices(self.all_directions, self.w)[0]
            #For loop to specify how many fibers to be generated when it is called
            for i in range(1):
                # pass
                self.ECM_Gen(r.pixel.x, r.pixel.y, r.pixel.z, init_direction, Fiber_Length,self.NECM)
        # print('Cell numer', len(self.cell_list_by_type(self.AVEC)), len(list(pixel_list)))


    def ECM_Gen(self,x,y,z,direction,Fiber_Length,fiber_type):
        for i in range(Fiber_Length):
            # Length = 3 pixels
            # direction = random.randint(1, 8)  # When attempted for random direction

            if ((x + (self.direction_dict[direction][0]) > 0) and
                    (x + (self.direction_dict[direction][0]) < self.dim.x) and
                    (y + (self.direction_dict[direction][1]) > 0) and
                    (y + (self.direction_dict[direction][1]) < self.y_limit) and
                    (z + (self.direction_dict[direction][2] * i) > 0) and
                    (z + (self.direction_dict[direction][2] * i) < self.dim.z)):

                if i % 2 == 0:  # This parameter changes the direction of fiber curves after some
                    if not (direction == 30 or direction == 40):

                        #Creating left and right directions
                        left_direction = direction - 1
                        right_direction = direction + 1
                        if right_direction == 9 or right_direction == 19 or right_direction == 29:
                            right_direction -= 8
                        if left_direction == 0 or left_direction == 10 or left_direction == 20:
                            left_direction += 8

                        # Creating up and down directions
                        if direction < 10:
                            up_direction = direction + 10
                            down_direction = 30

                        elif 10 < direction < 20:
                            up_direction = direction + 10
                            down_direction = direction - 10

                        else: #20 < direction < 30
                            up_direction = 40
                            down_direction = direction - 10

                        direction = \
                            random.choices(
                                (left_direction, direction, right_direction, up_direction, down_direction),
                                (self.w[left_direction],
                                 self.w[direction],
                                 self.w[right_direction],
                                 self.w[up_direction],
                                 self.w[down_direction]))[0]

                    elif direction == 30:
                        direction = random.choices(
                            (direction, 1, 2, 3, 4, 5, 6, 7, 8),
                            (self.w[direction], self.w[1], self.w[2], self.w[3], self.w[4], self.w[5], self.w[6], self.w[7], self.w[8]))[0]

                    elif direction == 40:
                        direction = random.choices(
                            (direction, 21, 22, 23, 24, 25, 26, 27, 28),
                            (self.w[direction], self.w[21], self.w[22], self.w[23], self.w[24], self.w[25], self.w[26], self.w[27], self.w[28]))[0]

                x = x + (self.direction_dict[direction][0])
                y = y + (self.direction_dict[direction][1])
                z = z + (self.direction_dict[direction][2])
                cell = self.new_cell(fiber_type)
                self.cellField[x, y, z] = cell
                cell.targetVolume = cell.volume
                cell.lambdaVolume = 5000
