import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';
import { useMutation } from '@tanstack/react-query';
import { toast } from 'sonner';
import { runGWAS } from '@/lib/api';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Loader2 } from 'lucide-react';

interface GWASResults {
  plots: {
    manhattan_plot: string;
    qq_plot: string;
  },
  result_table_path: string;
}

const GwasPage = () => {
  const { t } = useTranslation();
  const [params, setParams] = useState({
    trait: 'trait_yield',
    n_pcs: 5,
  });
  const [gwasResults, setGwasResults] = useState<GWASResults | null>(null);

  const gwasMutation = useMutation({
    mutationFn: runGWAS,
    onSuccess: (data) => {
      toast.success("GWAS run completed successfully!");
      setGwasResults(data);
    },
    onError: (error) => {
       toast.error(t('common.error'), {
        description: error.message,
      });
    }
  });

  const handleParamChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setParams(prev => ({ ...prev, [name]: name === 'trait' ? value : parseInt(value) }));
  }

  const handleRunGwas = () => {
    // @ts-ignore
    gwasMutation.mutate(params);
  }

  return (
    <div className="space-y-6">
      <h1 className="text-3xl font-bold">{t('gwas.title')}</h1>

      <Card>
        <CardHeader>
          <CardTitle>GWAS Parameters</CardTitle>
          <CardDescription>{t('gwas.description')}</CardDescription>
        </CardHeader>
        <CardContent className="grid md:grid-cols-2 gap-6">
          <div className="space-y-2">
            <Label htmlFor="trait">Trait</Label>
            <Input id="trait" name="trait" type="text" value={params.trait} onChange={handleParamChange} />
          </div>
          <div className="space-y-2">
            <Label htmlFor="n_pcs">Number of Principal Components (PCs)</Label>
            <Input id="n_pcs" name="n_pcs" type="number" value={params.n_pcs} onChange={handleParamChange} />
          </div>
        </CardContent>
        <CardContent>
            <Button onClick={handleRunGwas} disabled={gwasMutation.isPending}>
                {gwasMutation.isPending && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
                {t('common.run')}
            </Button>
        </CardContent>
      </Card>

      {gwasMutation.isSuccess && gwasResults && (
        <Card>
          <CardHeader>
            <CardTitle>GWAS Results</CardTitle>
            <CardDescription>
              <Button variant="link" asChild>
                <a href={gwasResults.result_table_path} download>Download Results CSV</a>
              </Button>
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-8">
            <div>
              <h3 className="font-semibold mb-2">Manhattan Plot</h3>
              <img src={gwasResults.plots.manhattan_plot} alt="Manhattan Plot" className="border rounded-md w-full" />
            </div>
            <div>
              <h3 className="font-semibold mb-2">Q-Q Plot</h3>
              <img src={gwasResults.plots.qq_plot} alt="Q-Q Plot" className="border rounded-md max-w-md" />
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
};

export default GwasPage;
