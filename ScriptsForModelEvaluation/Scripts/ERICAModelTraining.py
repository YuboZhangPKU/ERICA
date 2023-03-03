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
import random

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
    print('Shape of data:', datas.shape)
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
    print('Shape of data:', datas.shape)
    print('%s\t%s\tData Preprocessing done > > > > > > > > > > \n' % ((time.asctime(time.localtime(time.time()))), src_path))
    return datas

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
                        help="Input multiple sequence alignment (MSA) file or file list, split by comma", required=True)
    parser.add_argument("-l", "--Label",
                        help="Input label file or file list, split by comma", required=True)
    parser.add_argument("-p", "--Population",
                        help="Number of populations involved in the analyses", choices=("4", "5"), required=True)
    parser.add_argument("-o", "--Output",
                        help="Output model name (default: a random ID)", required=False)
    parser.add_argument("-m", "--Model",
                        help="Using pre-trained moedl for fine-tuning", required=False)
    parser.add_argument("-e", "--Epochs",
                        help="Number of epochs", default='3', required=False)
    parser.add_argument("-b", "--Batch",
                        help="Batch size", default='12', required=False)
    parser.add_argument("-r", "--Rate",
                        help="Learning rate", default='0.0001', required=False)
    parser.add_argument("--Iteration",
                        help="Number of iterations", required=False)

    # load the parameters
    ARGS = parser.parse_args()
    MSAInput = ARGS.Input.split(',')
    LabelInput = ARGS.Label.split(',')
    PopulationCount = ARGS.Population
    if ARGS.Output:
        ResOutput = ARGS.Output
    else:
        if PopulationCount == "4":
            ResOutput = "four_taxon_model_" + ''.join([str(random.randint(0,9)) for i in range(4)])
        elif PopulationCount == "5":
            ResOutput = "five_taxon_model_" + ''.join([str(random.randint(0,9)) for i in range(4)])
    #LogFile = open(ResOutput+'.log', 'w')

    if ARGS.Model:
        ModelRestore = ARGS.Model
    n_epochs = int(ARGS.Epochs)
    batch_size = int(ARGS.Batch)
    learning_rate = float(ARGS.Rate)
    if ARGS.Iteration:   
        n_iterations = int(ARGS.Iteration)
    saving = 1
    training_ratio = 0.9

    if PopulationCount == "4":
        data = FourTaxonDataProcess(MSAInput[0])
    elif PopulationCount == "5":
        data = FiveTaxonDataProcess(MSAInput[0])
    if len(MSAInput) > 1:
        for MSAFile in MSAInput[1:]:
            if PopulationCount == "4":
                Tempdata = FourTaxonDataProcess(MSAFile)
            elif PopulationCount == "5":
                Tempdata = FiveTaxonDataProcess(MSAFile)
            if Tempdata.shape[1:] == data.shape[1:]:   
                data = np.concatenate((data, Tempdata))
            else:
                print('Inconsistent data. Program exit.')
                sys.exit()
        print('Shape of data:', data.shape)
        del Tempdata

    with open (LabelInput[0], 'r') as RawLabel:
        Labels = RawLabel.readlines()
        Labels = np.array([x.strip().split() for x in Labels], dtype=np.float)
    if len(LabelInput) > 1:
        for LabelFile in LabelInput[1:]:
            with open (LabelFile, 'r') as RawLabel:
                Label = RawLabel.readlines()
                Label = np.array([x.strip().split() for x in Label], dtype=np.float)
            if Label.shape[1] == Labels.shape[1]:   
                Labels = np.concatenate((Labels, Label))
            else:
                print('Inconsistent labels. Program exit.')
                sys.exit()
    print('Shape of labels:', Labels.shape)

    n_training = int(training_ratio * data.shape[0])
    random_index = np.arange(data.shape[0])
    np.random.shuffle(random_index)
    training_data = data[random_index][0:n_training]
    training_label = Labels[random_index][0:n_training]
    n_step = training_data.shape[0] // batch_size
    if not ARGS.Iteration:   
        n_iterations = n_step * n_epochs + 1  

    test_data = data[random_index][n_training:]
    test_label = Labels[random_index][n_training:]

    print('Shape of training dataset:', training_data.shape)
    print('Shape of test dataset:', test_data.shape)

    if test_data.shape[0] < batch_size:
        print('Too little test data. Program exit.')
        sys.exit()

    del data
    # load the Network
    if PopulationCount == "4":
        g = FourTaxonNetworks(batch_size=batch_size)
    elif PopulationCount == "5":
        g = FiveTaxonNetworks(batch_size=batch_size)    

    # model training
    with tf.compat.v1.Session(graph=g.graph) as sess:
        print('%s\tStart model training.' % (time.asctime(time.localtime(time.time()))))
        print("Iteration\tTime")
        #LogFile.write("Iteration\tTrianMAE\tTestMAE\n")
        global_step = tf.Variable(0, trainable=False)

        #loss = tf.compat.v1.keras.losses.MeanAbsoluteError(reduction=tf.keras.losses.Reduction.NONE)
        #mae = loss(g.y, g.logits)
        #mae = tf.reduce_mean(mae)
        mae = tf.compat.v1.losses.absolute_difference(g.y, g.logits)      
        tf.compat.v1.summary.scalar('MAE', mae)
        mse = tf.compat.v1.losses.mean_pairwise_squared_error(g.y, g.logits)
        tf.compat.v1.summary.scalar('MSE', mse)
        maeaddmse = add = tf.add(mae, mse)
        crossentropy = tf.keras.losses.categorical_crossentropy(g.y, g.logits)
        crossentropy = tf.reduce_mean(crossentropy)
        tf.compat.v1.summary.scalar('CrossEntropy', crossentropy)
        acc = tf.reduce_sum(tf.cast(tf.equal(tf.argmax(g.y, axis=-1),tf.argmax(g.logits, axis=-1)), tf.float32))/batch_size
        tf.compat.v1.summary.scalar('Accuracy', acc)
        mergeall = tf.compat.v1.summary.merge_all()
        train_summary_writer = tf.compat.v1.summary.FileWriter(ResOutput + '_summary/Train', graph=tf.compat.v1.get_default_graph())
        test_summary_writer = tf.compat.v1.summary.FileWriter(ResOutput + '_summary/Validate', graph=tf.compat.v1.get_default_graph())

        # using Mean Absolute Error or Mean Squared Error as the loss function
        #tf.compat.v1.disable_eager_execution()
        #optimizer = tf.compat.v1.train.GradientDescentOptimizer(learning_rate).minimize(mae)
        if PopulationCount == "4":
            optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate).minimize(maeaddmse, global_step=global_step)
        elif PopulationCount == "5":
            optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate).minimize(mae, global_step=global_step)
        
        init = tf.compat.v1.global_variables_initializer()
        sess.run(init)

        saver = tf.compat.v1.train.Saver() 
        if ARGS.Model:
            saver.restore(sess, ModelRestore)

        Flag = True
        Iteration = 0
        for epochs in range(n_epochs):
            if Flag:
                random_index = np.arange(training_data.shape[0])
                np.random.shuffle(random_index)
                random_data = training_data[random_index]
                random_label = training_label[random_index]
                for step in range(n_step):
                    xs = random_data[step * batch_size: (step+1) * batch_size]
                    ys = random_label[step * batch_size: (step+1) * batch_size]
                    _, summary = sess.run([optimizer, mergeall], feed_dict={g.x: xs, g.y: ys})
                    train_summary_writer.add_summary(summary, Iteration)

                    if Iteration % saving == 0:
                        #TrianMAE = sess.run(mae, feed_dict={g.x: xs, g.y: ys})
                        Test_index = np.random.randint(0, (test_data.shape[0]), batch_size)
                        xs = test_data[Test_index]
                        ys = test_label[Test_index]
                        summary2 = sess.run(mergeall, feed_dict={g.x: xs, g.y: ys})
                        test_summary_writer.add_summary(summary2, Iteration)
                        print("{}\t{}".format(Iteration, (time.asctime(time.localtime(time.time())))))
                        #TestMAE = sess.run(mae, feed_dict={g.x: xs, g.y: ys})
                        #print("{}\t{:.5f}\t{:.5f}\t{}".format(Iteration, TrianMAE, TestMAE, (time.asctime(time.localtime(time.time())))))
                        #LogFile.write("{}\t{:.5f}\t{:.5f}\n".format(Iteration, TrianMAE, TestMAE))
                    if Iteration >= n_iterations:
                        Flag = False
                        break
                    Iteration += 1
            else:
                break
        saver.save(sess, ResOutput)
        #LogFile.close()
        print('%s\tDone.' % (time.asctime(time.localtime(time.time()))))

