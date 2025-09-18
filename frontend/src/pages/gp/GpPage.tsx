import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';
import { useMutation } from '@tanstack/react-query';
import { toast } from 'sonner';
import { trainGP } from '@/lib/api';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Loader2 } from 'lucide-react';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table"

interface CVResults {
  r2: number;
  rmse: number;
  pearson_r: number;
}

interface GEBV {
  sample_id: string;
  gebv: number;
}

interface GPResults {
  run_id: number;
  model_id: string;
  cv_results: { [trait: string]: CVResults };
  gebvs: GEBV[];
}

const GpPage = () => {
  const { t } = useTranslation();
  const [trait, setTrait] = useState('trait_yield');
  const [gpResults, setGpResults] = useState<GPResults | null>(null);

  const gpMutation = useMutation({
    mutationFn: trainGP,
    onSuccess: (data) => {
      toast.success("GP model training completed successfully!");
      setGpResults(data);
    },
    onError: (error: any) => {
       toast.error(t('common.error'), {
        description: error.response?.data?.detail || error.message,
      });
    }
  });

  const handleRunTraining = () => {
    // In a real app, you would select train/candidate IDs.
    // For this MVP, we pass empty arrays and let the backend use all available data.
    gpMutation.mutate({
      traits: [trait],
      train_ids: [],
      candidate_ids: [],
      k_folds: 5
    });
  }

  const currentTraitCvResults = gpResults?.cv_results?.[trait];

  return (
    <div className="space-y-6">
      <h1 className="text-3xl font-bold">{t('gp.title')}</h1>

      <Card>
        <CardHeader>
          <CardTitle>GP Model Training</CardTitle>
          <CardDescription>{t('gp.description')}</CardDescription>
        </CardHeader>
        <CardContent className="grid md:grid-cols-3 gap-6 items-end">
          <div className="space-y-2">
            <Label htmlFor="trait">Trait to Train On</Label>
            <Input id="trait" name="trait" type="text" value={trait} onChange={(e) => setTrait(e.target.value)} />
          </div>
          <div>
            <Button onClick={handleRunTraining} disabled={gpMutation.isPending}>
                {gpMutation.isPending && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
                Train Model
            </Button>
          </div>
        </CardContent>
      </Card>

      {gpMutation.isSuccess && gpResults && currentTraitCvResults && (
        <div className="grid md:grid-cols-2 gap-6">
            <Card>
                <CardHeader>
                    <CardTitle>Cross-Validation Results</CardTitle>
                    <CardDescription>5-fold CV performance for trait: {trait}</CardDescription>
                </CardHeader>
                <CardContent className="grid grid-cols-3 gap-4 text-center">
                    <div><p className="text-sm text-muted-foreground">R-squared (R²)</p><p className="text-2xl font-bold">{currentTraitCvResults.r2.toFixed(3)}</p></div>
                    <div><p className="text-sm text-muted-foreground">Pearson's r</p><p className="text-2xl font-bold">{currentTraitCvResults.pearson_r.toFixed(3)}</p></div>
                    <div><p className="text-sm text-muted-foreground">RMSE</p><p className="text-2xl font-bold">{currentTraitCvResults.rmse.toFixed(3)}</p></div>
                </CardContent>
            </Card>
            <Card>
                <CardHeader>
                    <CardTitle>GEBV Predictions</CardTitle>
                    <CardDescription>Predictions for all samples using the trained model.</CardDescription>
                </CardHeader>
                <CardContent className="max-h-96 overflow-y-auto">
                    <Table>
                        <TableHeader>
                            <TableRow>
                                <TableHead>Sample ID</TableHead>
                                <TableHead className="text-right">GEBV</TableHead>
                            </TableRow>
                        </TableHeader>
                        <TableBody>
                            {gpResults.gebvs?.map(g => (
                                <TableRow key={g.sample_id}>
                                    <TableCell>{g.sample_id}</TableCell>
                                    <TableCell className="text-right">{g.gebv.toFixed(4)}</TableCell>
                                </TableRow>
                            ))}
                        </TableBody>
                    </Table>
                </CardContent>
            </Card>
        </div>
      )}
    </div>
  );
};

export default GpPage;
