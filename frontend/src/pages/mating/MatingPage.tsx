import { useTranslation } from 'react-i18next';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/card';

const MatingPage = () => {
  const { t } = useTranslation();

  return (
    <div className="space-y-6">
      <h1 className="text-3xl font-bold">{t('mating.title')}</h1>
      <Card>
        <CardHeader>
          <CardTitle>{t('mating.title')}</CardTitle>
          <CardDescription>{t('mating.description')}</CardDescription>
        </CardHeader>
        <CardContent>
          <p className="text-muted-foreground">
            Feature under construction. This section will allow you to:
          </p>
          <ul className="list-disc list-inside text-muted-foreground mt-2">
            <li>Set weights for multiple traits.</li>
            <li>Adjust the balance between genetic gain (alpha) and diversity (beta).</li>
            <li>Run a Genetic Algorithm to find optimal mating pairs.</li>
            <li>View a table of the best potential crosses.</li>
            <li>Visualize the gain vs. diversity trade-off on a scatter plot.</li>
          </ul>
        </CardContent>
      </Card>
    </div>
  );
};

export default MatingPage;
