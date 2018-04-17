#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>

#include <algorithm>
#include<bits/stdc++.h> 
#include <list>

using namespace std;

double pixel_size   = 0.00416667;
double min_x_center = -90.40409499;
double min_y_center = 19.72078911;

// number of workplaces to look at when selecting a workplace
int workplace_neighborhood = 1000;

// ratio for Mexico, according to World Bank
int student_teacher_ratio = 28;

// IPUMS values for EMPSTATD (detailed employment status) variable
vector<int> work_codes   = {110, 112, 113, 116, 120};
vector<int> home_codes   = {114, 200, 310, 320, 340, 390, 999};
vector<int> school_codes = {111, 330};
int child_code   = 0; // Used, inexplicably, for children 0-11
int school_age   = 5; // children under this stay home

// This is an optimization to cut down RAM by >50% and maintain
// speed by using lists with an index lookup rather than a dictionary
// to represent each person
//field_idx = dict(zip('hid age sex x y workid hh_serial pernum day_loc'.split(), range(9)))

double deg_to_rad(double degree) {
    return M_PIl*degree/180;
}

double rad_to_deg(double radians) {
    return 180*radians/M_PIl;
}


double haversine(double lon1, double lat1, double lon2, double lat2) {
    // Calculate the great circle distance between two points 
    // on the earth (specified in decimal degrees)

    // convert decimal degrees to radians 
    lon1 = deg_to_rad(lon1);
    lon2 = deg_to_rad(lon2);
    lat1 = deg_to_rad(lat1);
    lat2 = deg_to_rad(lat2);
    // haversine formula 
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;
    double a = pow(sin(dlat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dlon/2),2);
    double c = 2 * asin(sqrt(a));
    double km = 6371 * c;
    return km;
}

char lookup_location_code( int age,int  code) {
    char loc =  ' ';
   
// need to sort to use binary, (?)
sort(school_codes.begin(), school_codes.end());
sort(work_codes.begin(), work_codes.end() );
sort(home_codes.begin(), home_codes.end() );

    if (binary_search(school_codes.begin(),school_codes.end(), code) != school_codes.end() ) or (code == child_code and age >= school_age) {
    	loc = 's';

    }	else if (binary_search(work_codes.begin(), work_codes.end(), code) != work_codes.end() ) {
    	loc = 'w';

    }	else if (binary_search(home_codes.begin(), home_codes.end(), code) != home_codes.end() )  or (code == child_code and age < school_age) {
        loc = 'h';

    } else {
        cout <<  'Error:: encountered bad employment (EMPSTATD) code:' << "\n", code;
        exit();  

    }  return loc

}

//not sure if variable types are correct
//is this confusing since we are already using binarysearch above?

sort(val_list.begin(), val_list.end());

int binary_search(int val_list.begin(), int val_list.end() ,int val) {

//    print "val_list length:", len(val_list)
//    print "val, val_list[0][label], val_list[-1][label]:", val, val_list[0][label], val_list[-1][label]
 
	int val= 

    	if (val_list[0] >  val) {
        return 0;
	}  if (val_list[-1] < val) {
        return val_list.size() ;
	}

    int imax = val_list.size() ;
    imin = 0;

    while (imin < imax) {
        int imid = imin + (imax-imin)/2;
        if (val_list[imid] < val) {
            int imin = imid+1;
	}  else {
            imax = imid;
       	}
    }
   
    if ( (imax == imin) and val_list[imin] == val) {
        return imin;

//need extra help here

    }  else {
        if bound == 'lower'{
            return imin;
	}  if bound == 'upper' {
	} 	return imin + 1
	
}

double get_nearby_places(pxi, pyi, loc_type, workplaces_and_schools, num_loc_needed) {

    //print "searching for schools:", loc_type
  
    commute_range = -1;
    positions_found = 0;
    nearby_places = {}  // by index in workplaces_and_schools

//why are we not including position type at top 
    char pos_type= ' ';
 
    if loc_type == 'w'{
    pos_type = 'workers';

    } else {
        pos_type = 'students';
      } while (len(nearby_places) < num_loc_needed and len(nearby_places) < len(workplaces_and_schools)) or (positions_found <= 0) {
        nearby_places = {};
        commute_range += 1;
      }
        //print "\n\nEnvelope size:", 2*commute_range + 1, 'x', 2*commute_range + 1
	
        for x_val in range(pxi-commute_range, pxi+commute_range+1) {
            pos_xmin = binary_search(workplaces_and_schools, x_val, dict_val_lt, dict_val_gt, dict_val_eq, 'xi', 'lower');
            pos_xmax = binary_search(workplaces_and_schools, x_val+1, dict_val_lt, dict_val_gt, dict_val_eq, 'xi', 'upper');
	}    if pos_xmin == pos_xmax {
                continue;
	     }

            pos_ymin = double binary_search(workplaces_and_schools[pos_xmin:pos_xmax], pyi-commute_range, dict_val_lt, dict_val_gt, dict_val_eq, 'yi', 'lower');
            pos_ymax = double binary_search(workplaces_and_schools[pos_xmin:pos_xmax], pyi+commute_range+1, dict_val_lt, dict_val_gt, dict_val_eq, 'yi', 'upper');
            
            for i,w in enumerate(workplaces_and_schools[pos_xmin:pos_xmax][pos_ymin:pos_ymax]) {
                //raw_idx = pos_xmin + pos_ymin + i
	    }    if pos_type == 'students' {
                    positions_found += 1;
                    nearby_places.append(w['workid']);
                 }  else if w[pos_type] > 0 {
                        positions_found += w[pos_type];
                        nearby_places.append(w['workid']);
                    }
       
	   //print "local positions:", positions_found            
           //print "workplaces found:", len(nearby_places) 
           //print "Envelope size:", 2*commute_range + 1, 'x', 2*commute_range + 1

    return nearby_places, positions_found;

}


double select_nearest_school(px, py, nearby_places, workplaces_and_schools, workplace_lookup) {

    //nearby_places is a list of indeces for workplaces_and_schools
    
    closest_school, min_dist = double workplaces_and_schools[workplace_lookup[nearby_places[0]]], () 

   //() evaluates to greater than any number

    for workid in nearby_places {
        s = workplaces_and_schools[workplace_lookup[workid]];

        //print 'candidate:', workid, s
	
        d = double haversine(px, py, s['x'], s['y']);
    }    if d < min_dist {
            min_dist = d;
            closest_school = s;
    }   
    return closest_school, min_dist;
}

int x_to_col_num(x) {
    return int(round((x - min_x_center)/pixel_size))
}

int y_to_row_num(y) {
    return int(round((y - min_y_center)/pixel_size))
}

int  xy_cmp(a,b) {
    if a['xi'] < b['xi'] {
    }    return -1;
    else if a['xi'] > b['xi'] {
    }    return 1;
    else {
    }    if a['yi'] < b['yi'] {
    }   return -1;
	else if a['yi'] > b['yi'] {
	}    return 1;
        else {
	}    return 0;
}

int dict_val_lt(A,b,label) {
    if A[label] < b {
    }    return True;
    else {
    }    return False;
}

int dict_val_gt(A,b,label) {
    if A[label] > b {:
}   return True
    else {
    }    return False;
}

int dict_val_eq(A,b,label) {
    if A[label] == b {
    }    return True;
    else {
    }    return False;
}

int import_households(filename, hh_loc) {
    header = True;
    print "reading household locations";

    // At this point, only includes houses
    for line in file(filename) {

        ///
        id type x y x_ctr y_ctr
        1 house -89.6910220003 20.7015302009 -89.69159442 20.69995656
        2 house -89.700145483 20.6641877526 -89.69992776 20.66245653
        3 house -89.758249035 20.6954360352 -89.75826114 20.69578989
        4 house -89.6974011142 20.6551249287 -89.69576109 20.65412319
        ///

    }
        if header {
            header = False;
	}    continue

        p = line.split();
        if p[1] != 'house' {
            print 'Expecting only houses in the location file at this point.  Found:', line
            exit();
      int hh_loc[int(p[0])] = {'x':float(p[2]), 'y':float(p[3])};
	}  return;
}

double import_workplaces_and_schools(filename, workplaces_and_schools, current_max_loc_id) {
    header = True;
    print "reading workplaces & schools";
   int loc_id = current_max_loc_id + 1;
    total_raw_workers = 0;

    for line in file(filename) {
        
        ///
        W 1 -89.6264173747 20.9599660422
        W 2 -89.6116996054 20.964832419
        W 1 -89.6428676041 20.9771890368
        W 2 -89.6405575542 20.9678584161
        W 1 -89.6255278043 20.9746128005
        ///
    }
        if header {
            header = False;
	}    continue;

        p = line.split();
        w = {'workid':loc_id, 'type':p[0].lower(), 'x':float(p[2]), 'y':float(p[3]), 'workers':0, 'students':0};

        if w['type'] == 'w' {
            w['raw_workers'] = int(p[1])
            total_raw_workers += w['raw_workers']
	}  else {
            w['raw_workers'] = 0;
            
        w['xi'] = x_to_col_num(w['x']);
        w['yi'] = y_to_row_num(w['y']);
        workplaces_and_schools.append(w);
        loc_id += 1;
        
    workplaces_and_schools.sort(cmp=xy_cmp);

    workplace_lookup = dict();
	}
    for i in range(len(workplaces_and_schools)) {
        workplace_lookup[workplaces_and_schools[i]['workid']] = i;
    }    return total_raw_workers, workplace_lookup;
}



int import_population(filename, pop, pop_ids, day_loc_ctr) {
    /// This function imports a preliminary version of the population file
    that doesn't have people assigned to workplaces yet (perhaps among
    other things.)///

    print "reading population";
    header = True;
    i = -1;

    for line in file(filename) {
        i += 1
        //if i % 10000 == 0:
        //    print i
        
	///
        pid hid age sex gridx gridy workid hh_serial pernum empstatd
        1 1 31 1 0 0 -1 2748179000 1 110
        2 1 29 2 0 0 -1 2748179000 2 110
        3 1 10 2 0 0 -1 2748179000 3 0
        4 2 32 1 0 0 -1 2748114000 1 110
        ///
	
        if header {
//            fo.write(line)
            header = False;
	}    continue;

        p =int  map(int, line.split())

        hid = p[1];
        age = p[2];
        day_loc = double lookup_location_code(age, p[9]);
        day_loc_ctr[day_loc] += 1;
        //           hid, age, sex, x, y, workid, hh_serial, pernum, day_loc
        pop[p[0]] = double p[1:4] + [hh_loc[hid]['x'], hh_loc[hid]['y'], -1, p[7], p[8], day_loc]
        //pop[p[0]] = {'hid':hid, 'age':age, 'sex':p[3], 'x':hh_loc[hid]['x'], 'y':hh_loc[hid]['y'], 'day_loc':day_loc, 'workid':-1}
	
        pop_ids.append(p[0]);

    print "Population size:", len(pop);
    print;
    print "Total number of workers (IPUMS):", day_loc_ctr['w'];
    print "Total number of students (IPUMS):", day_loc_ctr['s'];
    print "Total number of homebodies (IPUMS):", day_loc_ctr['h'];
    print;
    return;
    }

def send_kids_to_school(pop, pop_ids, workplaces_and_schools, total_raw_workers, workplace_lookup):
    schools = [location for location in workplaces_and_schools if location['type'] == 's']

    students_allocated = 0
    for pid in pop_ids:
        loc_type = pop[pid][field_idx['day_loc']]
        if loc_type != 's':
            continue
        #print "looking at pid: ", pid
        num_loc_needed = 1

        px, py = pop[pid][field_idx['x']], pop[pid][field_idx['y']]
        pxi, pyi = x_to_col_num(px), y_to_row_num(py)

        nearby_places, positions_found = get_nearby_places(pxi, pyi, loc_type, schools, num_loc_needed)
        #for s in nearby_places:
        #    print s

        school, distance = select_nearest_school(px, py, nearby_places, workplaces_and_schools, workplace_lookup)
        # 'school' is an element in the workplaces_and_schools list
        pop[pid][field_idx['workid']] = school['workid'] 
        school['students'] += 1
        school['raw_workers'] += 1.0/student_teacher_ratio
        total_raw_workers     += 1.0/student_teacher_ratio
        students_allocated += 1
        if students_allocated % 10000 == 0:
            print "students sent to school:", students_allocated
        #print px, py, school, distance
    return total_raw_workers

def choose_workplace(px, py, nearby_places, workplaces_and_schools, workplace_lookup):
    raw_weights = [0.0 for i in range(len(nearby_places))]
    for i, workid in enumerate(nearby_places):
        w = workplaces_and_schools[workplace_lookup[workid]]
        dist = haversine(px, py, w['x'], w['y'])
        size = w['workers']
        raw_weights[i] = size / dist**2

    # normalize weights
    probs = []
    total = sum(raw_weights)
    for wt in raw_weights:
        probs.append(wt/total)
    
    r = random()
    for i,p in enumerate(probs):
        if r < p:
            return workplaces_and_schools[workplace_lookup[nearby_places[i]]]
        r -= p

    return workplaces_and_schools[workplace_lookup[nearby_places[-1]]]


# Import household location data
hh_loc = dict()
import_households('locations-yucatan.txt', hh_loc)

# Import workplace and school location data
#
# We are using workplace size (# employees) as a weight, but ignoring
# school size data currently, as we feel it is more realistic to send
# students to the nearest school.
workplaces_and_schools = []
max_loc_id = max(hh_loc.keys())
total_raw_workers, workplace_lookup = import_workplaces_and_schools('schools_and_workplaces.out', workplaces_and_schools, max_loc_id)

# Import population data
pop = OrderedDict()
pop_ids = []
day_loc_ctr = {'h':0, 'w':0, 's':0} # home / work / school
import_population('population-yucatan.txt', pop, pop_ids, day_loc_ctr)

print "Total number of non-teacher jobs (DENUE):", total_raw_workers
print "Student:Teacher ratio (World Bank):", student_teacher_ratio
print "Total number of needed teachers:", day_loc_ctr['s']/student_teacher_ratio
print "Total raw number of jobs (DENUE + needed teachers):", total_raw_workers + day_loc_ctr['s']/student_teacher_ratio

shuffle(pop_ids)

# Send kids to nearest school
# Needs to happen first so we know how many teachers are needed
total_raw_workers = send_kids_to_school(pop, pop_ids, workplaces_and_schools, total_raw_workers, workplace_lookup)

# normalize workplace sizes
employment_rescaling_factor = day_loc_ctr['w'] / float(total_raw_workers)
for place in workplaces_and_schools:
    place['workers'] = place['raw_workers'] * employment_rescaling_factor 

# Filehandle for file we're going to write
fo = file('population-yucatan_no_copy.txt','w')
fo.write('pid hid age sex hh_serial pernum workid\n')

# Make a copy so we can delete places from the original data structure as they fill up
W_AND_S_COPY = deepcopy(workplaces_and_schools)
ctr = 0 

# For each person
for pid in pop.keys():
    person = pop[pid]
    loc_type = person[field_idx['day_loc']]
    # Have them stay at home if they don't work or go to school
    if loc_type == 'h':
        person[field_idx['workid']] = person[field_idx['hid']] # person stays home
    # If they work, probabilistically choose a workplace based on where they live
    # and how many positions are available at each workplace
    elif loc_type == 'w':
        px, py = person[field_idx['x']], person[field_idx['y']]
        pxi, pyi = x_to_col_num(px), y_to_row_num(py)

        nearby_places, positions_found = get_nearby_places(pxi, pyi, loc_type, workplaces_and_schools, workplace_neighborhood)
        #print len(nearby_places), positions_found 
        workplace = choose_workplace(px, py, nearby_places, W_AND_S_COPY, workplace_lookup)
        workplace['workers'] -= 1 # remove one available job
        # If the selected workplace no longer has openings,
        # remove it from the list, so we don't have to consider it again
        if workplace['workers'] <= 0:
            for i,v in enumerate(workplaces_and_schools):
                if v['workid'] == workplace['workid']:
                    del workplaces_and_schools[i]
                    break
        person[field_idx['workid']] = workplace['workid'] # assign worker

    # Students already have the "workid" (prob should be called day_loc_id to avoid
    # confusion), so we don't have to do much for them, just output their info
    fo.write(' '.join(map(str,[
                       pid, 
                       person[field_idx['hid']],  
                       person[field_idx['age']],  
                       person[field_idx['sex']], 
                       person[field_idx['hh_serial']],
                       person[field_idx['pernum']],
                       person[field_idx['workid']],
                      ])) + '\n') 
    ctr += 1
    if ctr % 1000 == 0:
        print "placed", ctr, "people"
        print "workplace list size:", len(workplaces_and_schools)
    '''
    pid hid age sex hh_serial pernum workid 
    1 1 31 1 0 0 -1 2748179000 1 110
    2 1 29 2 0 0 -1 2748179000 2 110
    3 1 10 2 0 0 -1 2748179000 3 0
    4 2 32 1 0 0 -1 2748114000 1 110
    '''
    
// this is our test function 

int main() {
    cout << haversine(0, 0, -91, 21) << endl;
}///
