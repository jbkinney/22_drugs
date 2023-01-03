from post_process_input import *
import glob

# Concatenate report files for each split fasq file
report_files = glob.glob(os.path.join(report_file_dir, "*.txt"))
df       = (pd.read_csv(f, sep='\t') for f in report_files)
conc_df  = pd.concat(df, ignore_index=True)

# Make the report for each LID and sample
report_df = conc_df.groupby(['LID','sample']).agg(\
                {'successful reads':'sum',
                 'total number of lines':'sum'}).reset_index()
# rename the columns
report_df.columns = ['LID', 'sample', 'sample reads', 'LID reads']
# reorder the columns
report_df = report_df[['LID', 'LID reads', 'sample', 'sample reads']]
# make the percentage of success column
report_df['sample pct %'] = 100*report_df['sample reads']/report_df['LID reads']

# Change the formatting option to show just two floating points
# Change the formatting option to show thousand commas for int
report_df['sample pct %'] = report_df['sample pct %'].map(lambda x: '%2.2f' % x)
report_df['LID reads'] = report_df['LID reads'].map(lambda x: '{:,}'.format(x))
report_df['sample reads'] = report_df['sample reads'].map(lambda x: '{:,}'.format(x))


# Write the report in comma separated format for github
report_df.to_csv(report_file, sep=',', index=False)
