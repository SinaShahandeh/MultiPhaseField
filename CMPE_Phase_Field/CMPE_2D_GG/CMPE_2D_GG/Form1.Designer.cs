namespace CMPE_2D_GG
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.btnStart = new System.Windows.Forms.Button();
            this.label1 = new System.Windows.Forms.Label();
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            this.label2 = new System.Windows.Forms.Label();
            this.btnSave = new System.Windows.Forms.Button();
            this.folderBrowserDialog1 = new System.Windows.Forms.FolderBrowserDialog();
            this.lblSave = new System.Windows.Forms.Label();
            this.lbltn = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.pictureBox1 = new System.Windows.Forms.PictureBox();
            this.txtt0 = new System.Windows.Forms.TextBox();
            this.label4 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.txttimesteps = new System.Windows.Forms.TextBox();
            this.label6 = new System.Windows.Forms.Label();
            this.txtsaveres = new System.Windows.Forms.TextBox();
            this.label7 = new System.Windows.Forms.Label();
            this.txtSaveGsize = new System.Windows.Forms.TextBox();
            this.label8 = new System.Windows.Forms.Label();
            this.txtL = new System.Windows.Forms.TextBox();
            this.txtm = new System.Windows.Forms.TextBox();
            this.label9 = new System.Windows.Forms.Label();
            this.txtkappa = new System.Windows.Forms.TextBox();
            this.label10 = new System.Windows.Forms.Label();
            this.txtPz = new System.Windows.Forms.TextBox();
            this.label11 = new System.Windows.Forms.Label();
            this.txtmboxsize = new System.Windows.Forms.TextBox();
            this.label12 = new System.Windows.Forms.Label();
            this.txtnboxsize = new System.Windows.Forms.TextBox();
            this.label13 = new System.Windows.Forms.Label();
            this.txtDelt = new System.Windows.Forms.TextBox();
            this.label14 = new System.Windows.Forms.Label();
            this.txtDelx = new System.Windows.Forms.TextBox();
            this.label15 = new System.Windows.Forms.Label();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).BeginInit();
            this.SuspendLayout();
            // 
            // btnStart
            // 
            this.btnStart.BackgroundImageLayout = System.Windows.Forms.ImageLayout.Zoom;
            this.btnStart.Cursor = System.Windows.Forms.Cursors.Default;
            this.btnStart.FlatAppearance.BorderSize = 0;
            this.btnStart.Location = new System.Drawing.Point(132, 503);
            this.btnStart.Name = "btnStart";
            this.btnStart.Size = new System.Drawing.Size(131, 27);
            this.btnStart.TabIndex = 0;
            this.btnStart.Text = "Start Simulation";
            this.btnStart.UseVisualStyleBackColor = true;
            this.btnStart.Click += new System.EventHandler(this.btnStart_Click);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Font = new System.Drawing.Font("Microsoft Sans Serif", 10F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label1.Location = new System.Drawing.Point(27, 28);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(137, 17);
            this.label1.TabIndex = 1;
            this.label1.Text = "Input Parameters:";
            // 
            // openFileDialog1
            // 
            this.openFileDialog1.FileName = "openFileDialog1";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(27, 142);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(113, 13);
            this.label2.TabIndex = 2;
            this.label2.Text = "Save Result Directory:";
            // 
            // btnSave
            // 
            this.btnSave.Location = new System.Drawing.Point(146, 137);
            this.btnSave.Name = "btnSave";
            this.btnSave.Size = new System.Drawing.Size(75, 23);
            this.btnSave.TabIndex = 3;
            this.btnSave.Text = "Browse";
            this.btnSave.UseVisualStyleBackColor = true;
            this.btnSave.Click += new System.EventHandler(this.btnSave_Click);
            // 
            // lblSave
            // 
            this.lblSave.AutoSize = true;
            this.lblSave.Location = new System.Drawing.Point(230, 142);
            this.lblSave.Name = "lblSave";
            this.lblSave.Size = new System.Drawing.Size(91, 13);
            this.lblSave.TabIndex = 4;
            this.lblSave.Text = "[Saving Directory]";
            // 
            // lbltn
            // 
            this.lbltn.AutoSize = true;
            this.lbltn.Location = new System.Drawing.Point(355, 510);
            this.lbltn.Name = "lbltn";
            this.lbltn.Size = new System.Drawing.Size(13, 13);
            this.lbltn.TabIndex = 5;
            this.lbltn.Text = "0";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(286, 510);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(58, 13);
            this.label3.TabIndex = 6;
            this.label3.Text = "Time Step:";
            // 
            // pictureBox1
            // 
            this.pictureBox1.Location = new System.Drawing.Point(415, 175);
            this.pictureBox1.Name = "pictureBox1";
            this.pictureBox1.Size = new System.Drawing.Size(416, 381);
            this.pictureBox1.TabIndex = 7;
            this.pictureBox1.TabStop = false;
            // 
            // txtt0
            // 
            this.txtt0.Location = new System.Drawing.Point(301, 175);
            this.txtt0.Name = "txtt0";
            this.txtt0.Size = new System.Drawing.Size(67, 20);
            this.txtt0.TabIndex = 8;
            this.txtt0.Text = "10";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(27, 179);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(175, 13);
            this.label4.TabIndex = 9;
            this.label4.Text = "Initial time steps before optimization ";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(27, 214);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(219, 13);
            this.label5.TabIndex = 11;
            this.label5.Text = "Total number of time steps to finish simulation";
            // 
            // txttimesteps
            // 
            this.txttimesteps.Location = new System.Drawing.Point(301, 210);
            this.txttimesteps.Name = "txttimesteps";
            this.txttimesteps.Size = new System.Drawing.Size(67, 20);
            this.txttimesteps.TabIndex = 10;
            this.txttimesteps.Text = "20000";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(27, 249);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(267, 13);
            this.label6.TabIndex = 13;
            this.label6.Text = "Time step intervals for  saving the microstructure results";
            // 
            // txtsaveres
            // 
            this.txtsaveres.Location = new System.Drawing.Point(301, 245);
            this.txtsaveres.Name = "txtsaveres";
            this.txtsaveres.Size = new System.Drawing.Size(67, 20);
            this.txtsaveres.TabIndex = 12;
            this.txtsaveres.Text = "500";
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(27, 284);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(245, 13);
            this.label7.TabIndex = 15;
            this.label7.Text = "Time step intervals for  saving the grain size results";
            // 
            // txtSaveGsize
            // 
            this.txtSaveGsize.Location = new System.Drawing.Point(301, 280);
            this.txtSaveGsize.Name = "txtSaveGsize";
            this.txtSaveGsize.Size = new System.Drawing.Size(67, 20);
            this.txtSaveGsize.TabIndex = 14;
            this.txtSaveGsize.Text = "100";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(88, 390);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(13, 13);
            this.label8.TabIndex = 16;
            this.label8.Text = "L";
            // 
            // txtL
            // 
            this.txtL.Location = new System.Drawing.Point(107, 386);
            this.txtL.Name = "txtL";
            this.txtL.Size = new System.Drawing.Size(51, 20);
            this.txtL.TabIndex = 17;
            this.txtL.Text = "1";
            // 
            // txtm
            // 
            this.txtm.Location = new System.Drawing.Point(212, 385);
            this.txtm.Name = "txtm";
            this.txtm.Size = new System.Drawing.Size(51, 20);
            this.txtm.TabIndex = 19;
            this.txtm.Text = "2.0";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(192, 389);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(15, 13);
            this.label9.TabIndex = 18;
            this.label9.Text = "m";
            // 
            // txtkappa
            // 
            this.txtkappa.Location = new System.Drawing.Point(317, 385);
            this.txtkappa.Name = "txtkappa";
            this.txtkappa.Size = new System.Drawing.Size(51, 20);
            this.txtkappa.TabIndex = 21;
            this.txtkappa.Text = "3.0";
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Font = new System.Drawing.Font("Times New Roman", 12F, System.Drawing.FontStyle.Italic, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label10.Location = new System.Drawing.Point(293, 386);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(17, 19);
            this.label10.TabIndex = 20;
            this.label10.Text = "κ";
            // 
            // txtPz
            // 
            this.txtPz.Location = new System.Drawing.Point(317, 420);
            this.txtPz.Name = "txtPz";
            this.txtPz.Size = new System.Drawing.Size(51, 20);
            this.txtPz.TabIndex = 23;
            this.txtPz.Text = "0.0";
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(291, 428);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(19, 13);
            this.label11.TabIndex = 22;
            this.label11.Text = "Pz";
            this.label11.Click += new System.EventHandler(this.label11_Click);
            // 
            // txtmboxsize
            // 
            this.txtmboxsize.Location = new System.Drawing.Point(247, 314);
            this.txtmboxsize.Name = "txtmboxsize";
            this.txtmboxsize.Size = new System.Drawing.Size(51, 20);
            this.txtmboxsize.TabIndex = 25;
            this.txtmboxsize.Text = "500";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(27, 318);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(119, 13);
            this.label12.TabIndex = 24;
            this.label12.Text = "System size (grid points)";
            // 
            // txtnboxsize
            // 
            this.txtnboxsize.Location = new System.Drawing.Point(317, 314);
            this.txtnboxsize.Name = "txtnboxsize";
            this.txtnboxsize.Size = new System.Drawing.Size(51, 20);
            this.txtnboxsize.TabIndex = 26;
            this.txtnboxsize.Text = "500";
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.Location = new System.Drawing.Point(300, 318);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(14, 13);
            this.label13.TabIndex = 27;
            this.label13.Text = "X";
            // 
            // txtDelt
            // 
            this.txtDelt.Location = new System.Drawing.Point(317, 350);
            this.txtDelt.Name = "txtDelt";
            this.txtDelt.Size = new System.Drawing.Size(51, 20);
            this.txtDelt.TabIndex = 31;
            this.txtDelt.Text = "0.1";
            // 
            // label14
            // 
            this.label14.AutoSize = true;
            this.label14.Location = new System.Drawing.Point(290, 354);
            this.label14.Name = "label14";
            this.label14.Size = new System.Drawing.Size(20, 13);
            this.label14.TabIndex = 30;
            this.label14.Text = "Δ t";
            // 
            // txtDelx
            // 
            this.txtDelx.Location = new System.Drawing.Point(212, 350);
            this.txtDelx.Name = "txtDelx";
            this.txtDelx.Size = new System.Drawing.Size(51, 20);
            this.txtDelx.TabIndex = 29;
            this.txtDelx.Text = "1.0";
            // 
            // label15
            // 
            this.label15.AutoSize = true;
            this.label15.Location = new System.Drawing.Point(180, 354);
            this.label15.Name = "label15";
            this.label15.Size = new System.Drawing.Size(22, 13);
            this.label15.TabIndex = 28;
            this.label15.Text = "Δ x";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(843, 568);
            this.Controls.Add(this.txtDelt);
            this.Controls.Add(this.label14);
            this.Controls.Add(this.txtDelx);
            this.Controls.Add(this.label15);
            this.Controls.Add(this.label13);
            this.Controls.Add(this.txtnboxsize);
            this.Controls.Add(this.txtmboxsize);
            this.Controls.Add(this.label12);
            this.Controls.Add(this.txtPz);
            this.Controls.Add(this.label11);
            this.Controls.Add(this.txtkappa);
            this.Controls.Add(this.label10);
            this.Controls.Add(this.txtm);
            this.Controls.Add(this.label9);
            this.Controls.Add(this.txtL);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.txtSaveGsize);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.txtsaveres);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.txttimesteps);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.txtt0);
            this.Controls.Add(this.pictureBox1);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.lbltn);
            this.Controls.Add(this.lblSave);
            this.Controls.Add(this.btnSave);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.btnStart);
            this.Name = "Form1";
            this.Text = "CMPE Multi Phase Field Simulator - 2D Grain Growth";
            this.Load += new System.EventHandler(this.Form1_Load);
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button btnStart;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.OpenFileDialog openFileDialog1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Button btnSave;
        private System.Windows.Forms.FolderBrowserDialog folderBrowserDialog1;
        private System.Windows.Forms.Label lblSave;
        private System.Windows.Forms.Label lbltn;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.PictureBox pictureBox1;
        private System.Windows.Forms.TextBox txtt0;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.TextBox txttimesteps;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.TextBox txtsaveres;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.TextBox txtSaveGsize;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.TextBox txtL;
        private System.Windows.Forms.TextBox txtm;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.TextBox txtkappa;
        private System.Windows.Forms.Label label10;
        private System.Windows.Forms.TextBox txtPz;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.TextBox txtnboxsize;
        private System.Windows.Forms.Label label13;
        private System.Windows.Forms.TextBox txtDelt;
        private System.Windows.Forms.Label label14;
        private System.Windows.Forms.TextBox txtDelx;
        private System.Windows.Forms.Label label15;
        private System.Windows.Forms.TextBox txtmboxsize;
    }
}

