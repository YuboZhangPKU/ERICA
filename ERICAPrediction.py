# _*_ coding: UTF-8 _*_
# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: ZHU QINGJIE, QIU JIWEN, ZHANG YUBO

    Version-01:
        2020-01  Inferring evolutionary relationship from multiple sequence alignment for four population

    Version-02:
        2020-08  Inferring evolutionary relationship from multiple sequence alignment for four and five population       
    """
# Version information END ----------------------------------------------------


import argparse
from multiprocessing import Process
import numpy as np
import os
from pathlib import Path
import sys
import tensorflow as tf
import tensorflow.compat.v1.layers as L
from tensorflow.python.keras.layers import MaxPool2D
from tensorflow.python.keras.layers import Activation
from tensorflow.python.keras.layers import Concatenate
from tensorflow.python.keras.layers import Conv2D
import time


TOKEN = {'G': 0, 'T': 1, 'A': 2, 'C': 3, 'N': 4, '-': 5, 'g': 0, 't': 1, 'a': 2, 'c': 3, 'n': 4}
###############################################################################
# function part
###############################################################################
# data preprocessing, sequence encoded in the one-hot form
def FourTaxonCalculateDiff(datas):
    mask0 = np.full(datas.shape, 0, dtype = np.uint8)
    mask1 = np.full(datas.shape, 1, dtype = np.uint8)
    mask2 = np.full(datas.shape, 2, dtype = np.uint8)
    mask3 = np.full(datas.shape, 3, dtype = np.uint8)
    x0 = np.equal(datas, mask0).astype(np.uint8)
    x1 = np.equal(datas, mask1).astype(np.uint8)
    x2 = np.equal(datas, mask2).astype(np.uint8)
    x3 = np.equal(datas, mask3).astype(np.uint8)
    datas = np.concatenate([x0, x1, x2, x3], axis=0)
    return datas


def FiveTaxonCalculateDiff(datas):
    mask0 = np.full(datas.shape, 0, dtype = np.uint8)
    mask1 = np.full(datas.shape, 1, dtype = np.uint8)
    mask2 = np.full(datas.shape, 2, dtype = np.uint8)
    mask3 = np.full(datas.shape, 3, dtype = np.uint8)
    x0 = np.equal(datas, mask0).astype(np.uint8)
    x1 = np.equal(datas, mask1).astype(np.uint8)
    x2 = np.equal(datas, mask2).astype(np.uint8)
    x3 = np.equal(datas, mask3).astype(np.uint8)
    datas = np.stack([x0, x1, x2, x3], axis=1)
    datas = np.reshape(datas, (160,5000))
    return datas


def FourTaxonDataProcess(src_path):
    print('%s\t%s\tData Preprocessing > > > > > > > > > > \n' % ((time.asctime(time.localtime(time.time()))), src_path))
    lines = open(src_path, 'r').readlines()
    def FnToken(x):
        for key in TOKEN:
            x = x.replace(key, str(TOKEN[key]))
        return x
    lines = [FnToken(x.strip()) for x in lines]
    lines = np.array([list(x) for x in lines], dtype=np.uint8)
    print('Shape of MSA:', lines.shape)
    data_len = lines.shape[1]
    result = data_len%5000
    if result:
        datas = np.array([lines[:, i:i + 5000] for i in range(0, data_len, 5000)][:-1])
    else:
        datas = np.array([lines[:, i:i + 5000] for i in range(0, data_len, 5000)])
    datas = np.array([FourTaxonCalculateDiff(i) for i in datas], dtype=np.uint8)
    print('Shape of datas:', datas.shape)
    print('%s\t%s\tData Preprocessing done > > > > > > > > > > \n' % ((time.asctime(time.localtime(time.time()))), src_path))
    return datas


def FiveTaxonDataProcess(src_path):
    print('%s\t%s\tData Preprocessing > > > > > > > > > > \n' % ((time.asctime(time.localtime(time.time()))), src_path))
    lines = open(src_path, 'r').readlines()
    def FnToken(x):
        for key in TOKEN:
            x = x.replace(key, str(TOKEN[key]))
        return x
    lines = [FnToken(x.strip()) for x in lines]
    lines = np.array([list(x) for x in lines], dtype=np.uint8)
    print('Shape of MSA:', lines.shape)
    data_len = lines.shape[1]
    result = data_len%5000
    if result:
        datas = np.array([lines[:, i:i + 5000] for i in range(0, data_len, 5000)][:-1])
    else:
        datas = np.array([lines[:, i:i + 5000] for i in range(0, data_len, 5000)])
    datas = np.array([FiveTaxonCalculateDiff(i) for i in datas], dtype=np.uint8)
    print('Shape of datas:', datas.shape)
    print('%s\t%s\tData Preprocessing done > > > > > > > > > > \n' % ((time.asctime(time.localtime(time.time()))), src_path))
    return datas


# data prediction by trained networks
def FourTaxonDataPrediction(MSA, res, g, Model):
    with tf.compat.v1.Session(graph=g.graph) as sess:
        saver = tf.compat.v1.train.Saver()
        saver.restore(sess, Model)
        data = FourTaxonDataProcess(MSA)
        print('%s\tPredictions will be writen into %s.' % ((time.asctime(time.localtime(time.time()))), res))
        xs1 = [data[i:i + 1] for i in range(0, len(data), 1)]
        resf = open(res, 'w')
        for xs in xs1:
            prd = sess.run([g.logits], feed_dict={g.x: xs})
            for p in prd:
                p = p[0]
                resf.write('{:.4f}   {:.4f}  {:.4f}  \n'.format(p[0], p[1], p[2]))
        resf.close()
        print('%s\tDone.' % (time.asctime(time.localtime(time.time()))))


def FiveTaxonDataPrediction(MSA, res, g, Model):
    with tf.compat.v1.Session(graph=g.graph) as sess:
        saver = tf.compat.v1.train.Saver()
        saver.restore(sess, Model)
        data = FiveTaxonDataProcess(MSA)
        print('%s\tPredictions will be writen into %s.' % ((time.asctime(time.localtime(time.time()))), res))
        xs1 = [data[i:i + 1] for i in range(0, len(data), 1)]
        resf = open(res, 'w')
        for xs in xs1:
            prd = sess.run([g.logits], feed_dict={g.x: xs})
            for p in prd:
                p = p[0]
                resf.write('{:.4f}   {:.4f}  {:.4f}  {:.4f}   {:.4f}  {:.4f}   {:.4f}  {:.4f}  {:.4f}   {:.4f}  {:.4f}   {:.4f}  {:.4f}  {:.4f}   {:.4f}  \n'.format(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14]))
        resf.close()
        print('%s\tDone.' % (time.asctime(time.localtime(time.time()))))


# architectures of the nerual networks 
class FourTaxonNetworks():
    def conv_block(self, x, growth_rate):
        x1 = Activation('relu')(x)
        x1 = Conv2D(growth_rate, 1)(x1)
        x1 = Activation('relu')(x1)
        x1 = Conv2D(growth_rate, 3, padding='same')(x1)
        x = Concatenate(axis=3)([x, x1])
        return x   


    def dense_block(self, x, blocks):
        for i in range(blocks):
            x = self.conv_block(x, 5)
        return x


    def _residual_block_first(self, x, out_channel, strides=(2, 3)):
        in_channel = x.get_shape().as_list()[-1]
        if in_channel == out_channel:
            if strides == 1:
                shortcut = x
            else:
                shortcut = MaxPool2D(pool_size=strides, strides=strides, padding='same')(x)
        else:
            shortcut = Conv2D(out_channel, 1, strides=strides, padding='same')(x)
        x = Activation('relu')(x)
        x = Conv2D(out_channel, kernel_size=strides, strides=strides, padding='same')(x)
        x = Activation('relu')(x)
        x = Conv2D(out_channel, kernel_size=strides, padding='same')(x)
        x = Activation('relu')(x)
        x = Conv2D(out_channel, kernel_size=strides, padding='same')(x)
        x = x + shortcut
        return x


    def __init__(sf, tranning=False, predict_flag=False, batch_size=8):
        sf.graph = tf.Graph()
        with sf.graph.as_default():
            sf.x = tf.compat.v1.placeholder(tf.int32, [batch_size, 128, 5000])
            sf.y = tf.compat.v1.placeholder(tf.float32, [batch_size, 3])
            x = tf.cast(sf.x, tf.float32)
            y_one_hot = sf.y
            x = tf.expand_dims(x, -1)
            x = Conv2D(32, 1, padding='same')(x)
            blocks = 3
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(8, 1))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(4, 1))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(1, 2))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(1, 3))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(1, 3))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(1, 3))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(2, 3))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = Activation('relu')(x)
            x = tf.squeeze(x)
            x = tf.reshape(x, [batch_size, -1])
            x = L.dense(x, 512, activation=tf.nn.relu)
            sf.logits = L.dense(x, 3, activation=tf.nn.softmax)


class FiveTaxonNetworks():
    def conv_block(self, x, growth_rate):
        x1 = Activation('relu')(x)
        x1 = Conv2D(growth_rate, 1)(x1)
        x1 = Activation('relu')(x1)
        x1 = Conv2D(growth_rate, 3, padding='same')(x1)
        x = Concatenate(axis=3)([x, x1])
        return x   


    def dense_block(self, x, blocks):
        for i in range(blocks):
            x = self.conv_block(x, 8)
        return x


    def _residual_block_first(self, x, out_channel, strides=(2, 3)):
        in_channel = x.get_shape().as_list()[-1]
        if in_channel == out_channel:
            if strides == 1:
                shortcut = x
            else:
                shortcut = MaxPool2D(pool_size=strides, strides=strides, padding='same')(x)
        else:
            shortcut = Conv2D(out_channel, 1, strides=strides, padding='same')(x)
        x = Activation('relu')(x)
        x = Conv2D(out_channel, kernel_size=strides, strides=strides, padding='same')(x)
        x = Activation('relu')(x)
        x = Conv2D(out_channel, kernel_size=strides, padding='same')(x)
        x = Activation('relu')(x)
        x = Conv2D(out_channel, kernel_size=strides, padding='same')(x)
        x = x + shortcut
        return x


    def __init__(sf, tranning=False,predict_flag=False, batch_size=10):
        sf.graph = tf.Graph()
        with sf.graph.as_default():
            sf.x = tf.compat.v1.placeholder(tf.int32, [batch_size, 160, 5000])
            sf.y = tf.compat.v1.placeholder(tf.float32, [batch_size, 15])
            x = tf.cast(sf.x, tf.float32)
            y_one_hot = sf.y
            x = tf.expand_dims(x, -1)
            x = Conv2D(32, 1,  padding='same')(x)
            blocks = 3
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(4, 1))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel,strides=(8,1))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel,strides=(1,2))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel,strides=(1,3))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(1, 3))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(1, 3))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(2, 3))
            x = sf.dense_block(x, blocks)
            in_channel = x.get_shape().as_list()[-1]
            x = sf._residual_block_first(x, in_channel, strides=(2, 3))
            x = sf.dense_block(x, blocks)
            x = Activation('relu')(x)
            x = tf.squeeze(x)
            x = tf.reshape(x, [batch_size, -1])
            x = L.dense(x, 512, activation=tf.nn.relu)
            sf.logits = L.dense(x, 15, activation=tf.nn.softmax)


# read parameters
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evolutionary Relationship Inference using a CNN-based Approach")
    parser.add_argument("-i", "--Input",
                        help="Input multiple sequence alignment (MSA) file or directory", required=True)
    parser.add_argument("-o", "--Output",
                        help="Output file name", required=False)
    parser.add_argument("-p", "--Population",
                        help="Number of populations involved in the analyses", choices=("4", "5"), required=True)
    parser.add_argument("-t", "--Tasks", default="1",
                        help="Number of tasks running in parallel. Only used for multiple input files", required=False)
    parser.add_argument("-m", "--Model",
                        help="Model used for prediction", required=False)


    # load the parameters
    ARGS = parser.parse_args()
    MSAInput = ARGS.Input
    ResOutput = ARGS.Output
    PopulationCount = ARGS.Population


    # default models for four-taxon and five-taxon analyses 
    if PopulationCount == "4":
        Model = 'TrainedModels/four_taxon_model_319200'
        g = FourTaxonNetworks(predict_flag=True, batch_size=1)
    elif PopulationCount == "5":
        Model = 'TrainedModels/five_taxon_model_660600'    
        g = FiveTaxonNetworks(predict_flag=True, batch_size=1)
    if ARGS.Model:
        Model = ARGS.Model


    # input a directory
    if Path(MSAInput).is_dir():
        Tasks = int(ARGS.Tasks)
        MSAInputList = os.listdir(MSAInput)
        ResOutputList = []
        if ResOutput:
            if not os.path.exists(ResOutput):
                os.makedirs(ResOutput)
            for i in range(len(MSAInputList)):
                if MSAInputList[i][-4:] == '.txt':
                    ResOutputList.append(os.path.abspath(ResOutput) + '/' + MSAInputList[i][:-4] + '_res.txt')
                else:
                    ResOutputList.append(os.path.abspath(ResOutput) + '/' + MSAInputList[i] + '_res.txt')
                MSAInputList[i] = (os.path.abspath(MSAInput) + '/' + MSAInputList[i])
        else:
            for i in range(len(MSAInputList)):
                if MSAInputList[i][-4:] == '.txt':
                    ResOutputList.append(os.path.abspath(MSAInput) + '/' + MSAInputList[i][:-4] + '_res.txt')
                else:
                    ResOutputList.append(os.path.abspath(MSAInput) + '/' + MSAInputList[i] + '_res.txt')
                MSAInputList[i] = (os.path.abspath(MSAInput) + '/' + MSAInputList[i])
    

    # input an MSA file
    elif Path(MSAInput).is_file():
        Tasks = 1
        MSAInputList = [MSAInput]
        ResOutputList = []
        if ResOutput:
            ResOutputList = [ResOutput]
        else:
            if MSAInput[-4:] == '.txt':
                ResOutputList.append(os.path.abspath(MSAInput)[:-4] + '_res.txt')
            else:
                ResOutputList.append(os.path.abspath(MSAInput) + '_res.txt')
    else:
        print("There is no file or dir named %s. Program exit." % (MSAInput))
        sys.exit()


    # running prediction jobs in parallel
    print("-" * 30)
    index = 0
    flag = True
    while flag:
        processes = []
        for i in range(Tasks):
            try:
                if PopulationCount == "4":
                    p = Process(target=FourTaxonDataPrediction, args=(MSAInputList[index], ResOutputList[index], g, Model))
                elif PopulationCount == "5":
                    p = Process(target=FiveTaxonDataPrediction, args=(MSAInputList[index], ResOutputList[index], g, Model))                  
                p.start()
                processes.append(p)
                index += 1
            except:
                flag = False
                continue
        for p in processes:
            p.join()
