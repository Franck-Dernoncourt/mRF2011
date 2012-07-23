#! /usr/bin/env python
def build(bld):
	obj = bld.new_task_gen('cxx', 'program')
	obj.source = 'mrf.cpp'
	obj.includes = '. ../../'
	obj.uselib = 'BOOST BOOST_UNIT_TEST_FRAMEWORK EIGEN2 GSL BOOST_GRAPH BOOST_SYSTEM BOOST_THREAD BOOST_FILESYSTEM'
	obj.uselib_local = 'sferes2'
	obj.target = 'mrf'
	



