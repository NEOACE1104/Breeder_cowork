import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';
import { useMutation } from '@tanstack/react-query';
import { toast } from 'sonner';
import { uploadFile, runQC } from '@/lib/api';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Loader2 } from 'lucide-react';

interface QCResults {
  summary: {
    initial_samples: number;
    initial_snps: number;
    samples_after_qc: number;
    snps_after_qc: number;
    mean_maf: number;
  },
  plots: {
    maf_distribution: string;
  }
}

const DataQcPage = () => {
  const { t } = useTranslation();
  const [file, setFile] = useState<File | null>(null);
  const [qcParams, setQcParams] = useState({
    max_sample_missing: 0.2,
    max_snp_missing: 0.2,
    min_maf: 0.05,
    min_hwe_p: 1e-6,
  });
  const [qcResults, setQcResults] = useState<QCResults | null>(null);

  const uploadMutation = useMutation({
    mutationFn: uploadFile,
    onSuccess: (data) => {
      toast.success(t('data_qc.upload_success'), {
        description: `Loaded ${data.num_samples} samples and ${data.num_snps} SNPs.`,
      });
    },
    onError: (error) => {
      toast.error(t('common.error'), {
        description: error.message,
      });
    },
  });

  const qcMutation = useMutation({
    mutationFn: runQC,
    onSuccess: (data) => {
      toast.success("QC run completed successfully!");
      setQcResults(data);
    },
    onError: (error) => {
       toast.error(t('common.error'), {
        description: error.message,
      });
    }
  });

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files) {
      setFile(e.target.files[0]);
    }
  };

  const handleUpload = () => {
    if (file) {
      uploadMutation.mutate(file);
    }
  };

  const handleParamChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setQcParams(prev => ({...prev, [name]: parseFloat(value) }));
  }

  const handleRunQc = () => {
    qcMutation.mutate(qcParams);
  }

  return (
    <div className="space-y-6">
      <h1 className="text-3xl font-bold">{t('data_qc.title')}</h1>

      <Card>
        <CardHeader>
          <CardTitle>{t('data_qc.upload_title')}</CardTitle>
          <CardDescription>{t('data_qc.upload_desc')}</CardDescription>
        </CardHeader>
        <CardContent className="flex gap-4 items-center">
          <Input type="file" onChange={handleFileChange} className="max-w-xs" />
          <Button onClick={handleUpload} disabled={!file || uploadMutation.isPending}>
            {uploadMutation.isPending && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
            {t('data_qc.upload_button')}
          </Button>
        </CardContent>
      </Card>

      <Card>
        <CardHeader>
          <CardTitle>{t('data_qc.qc_params')}</CardTitle>
        </CardHeader>
        <CardContent className="grid md:grid-cols-2 gap-6">
          <div className="space-y-2">
            <Label htmlFor="max_sample_missing">Max Sample Missing Rate</Label>
            <Input id="max_sample_missing" name="max_sample_missing" type="number" value={qcParams.max_sample_missing} onChange={handleParamChange} step="0.05" />
          </div>
          <div className="space-y-2">
            <Label htmlFor="max_snp_missing">Max SNP Missing Rate</Label>
            <Input id="max_snp_missing" name="max_snp_missing" type="number" value={qcParams.max_snp_missing} onChange={handleParamChange} step="0.05" />
          </div>
          <div className="space-y-2">
            <Label htmlFor="min_maf">Min MAF</Label>
            <Input id="min_maf" name="min_maf" type="number" value={qcParams.min_maf} onChange={handleParamChange} step="0.01" />
          </div>
          <div className="space-y-2">
            <Label htmlFor="min_hwe_p">Min HWE P-value</Label>
            <Input id="min_hwe_p" name="min_hwe_p" type="number" value={qcParams.min_hwe_p} onChange={handleParamChange} step="1e-7" />
          </div>
        </CardContent>
        <CardContent>
            <Button onClick={handleRunQc} disabled={qcMutation.isPending || !uploadMutation.isSuccess}>
                {qcMutation.isPending && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
                {t('data_qc.run_qc')}
            </Button>
        </CardContent>
      </Card>

      {qcMutation.isSuccess && qcResults && (
        <Card>
          <CardHeader>
            <CardTitle>{t('data_qc.qc_results')}</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
             <div className="grid grid-cols-2 md:grid-cols-5 gap-4 text-center">
                <div><p className="text-sm text-muted-foreground">{t('data_qc.initial_samples')}</p><p className="text-2xl font-bold">{qcResults.summary.initial_samples}</p></div>
                <div><p className="text-sm text-muted-foreground">{t('data_qc.initial_snps')}</p><p className="text-2xl font-bold">{qcResults.summary.initial_snps}</p></div>
                <div><p className="text-sm text-muted-foreground">{t('data_qc.samples_after_qc')}</p><p className="text-2xl font-bold">{qcResults.summary.samples_after_qc}</p></div>
                <div><p className="text-sm text-muted-foreground">{t('data_qc.snps_after_qc')}</p><p className="text-2xl font-bold">{qcResults.summary.snps_after_qc}</p></div>
                <div><p className="text-sm text-muted-foreground">{t('data_qc.mean_maf')}</p><p className="text-2xl font-bold">{qcResults.summary.mean_maf.toFixed(3)}</p></div>
            </div>
            <div>
              <h3 className="font-semibold mb-2">MAF Distribution</h3>
              <img src={qcResults.plots.maf_distribution} alt="MAF Distribution Plot" className="border rounded-md" />
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
};

export default DataQcPage;
